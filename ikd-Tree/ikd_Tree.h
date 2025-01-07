#pragma once
#include <algorithm>
#include <chrono>
#include <math.h>
#include <memory>
#include <pcl/point_types.h>
#include <pthread.h>
#include <queue>
#include <stdio.h>
#include <time.h>
#include <unistd.h>

#define EPSS 1e-6
#define Minimal_Unbalanced_Tree_Size 10
#define Multi_Thread_Rebuild_Point_Num 1500
#define DOWNSAMPLE_SWITCH true
#define ForceRebuildPercentage 0.2
#define Q_LEN 1000000

using namespace std;

struct ikdTree_PointType {
  float x, y, z;
  ikdTree_PointType(float px = 0.0f, float py = 0.0f, float pz = 0.0f) {
    x = px;
    y = py;
    z = pz;
  }
};

struct BoxPointType {
  float vertex_min[3];
  float vertex_max[3];
};

enum operation_set {
  ADD_POINT,
  DELETE_POINT,
  DELETE_BOX,
  ADD_BOX,
  DOWNSAMPLE_DELETE,
  PUSH_DOWN
};

enum delete_point_storage_set {
  NOT_RECORD,
  DELETE_POINTS_REC,
  MULTI_THREAD_REC
};
// 自定义队列
template <typename T>
class MANUAL_Q {
 private:
  int head = 0, tail = 0, counter = 0;  // counter记录队列内元素个数
  T q[Q_LEN];                           // Q_LEN表示队列最大容量
  bool is_empty;

 public:
  void pop();
  T front();
  T back();
  void clear();
  void push(T op);
  bool empty();
  int size();
};

template <typename PointType>
class KD_TREE {
 public:
  using PointVector = vector<PointType, Eigen::aligned_allocator<PointType>>;
  using Ptr = shared_ptr<KD_TREE<PointType>>;
  struct KD_TREE_NODE {
    PointType point;
    uint8_t division_axis;  // 标记划分的坐标轴
    int TreeSize = 1;  // 以当前节点为根节点的subtree中所有节点的个数
    int invalid_point_num = 0;  // 失效点数量，即deleted标记为true的点个数
    int down_del_num = 0;                   // 下采样时删除点云个数
    bool point_deleted = false;             // 标记当前节点是否被删除
    bool tree_deleted = false;              // 整棵子树是否标记为删除
    bool point_downsample_deleted = false;  // 标记该节点下采样时被删除
    bool tree_downsample_deleted = false;  // 标记subtree下采样时被删除
    bool need_push_down_to_left =
        false;  // 左子树需要执行push_down操作（置true时）
    bool need_push_down_to_right =
        false;  // 右子树需要执行push_down操作（置false时）
    bool working_flag = false;  // 标记节点是否被占用，处于work状态
    float radius_sq;  // bbox的半径平方，最大最小点差值的一半的平方和
    pthread_mutex_t push_down_mutex_lock;  // 互斥锁
    float node_range_x[2], node_range_y[2],
        node_range_z
            [2];  // 以当前节点为根节点的subtree中所有节点，在三维空间中的最小bbox，对应的max_point、min_point
                  // // 剪枝加速都依赖该变量
    KD_TREE_NODE* left_son_ptr = nullptr;
    KD_TREE_NODE* right_son_ptr = nullptr;
    KD_TREE_NODE* father_ptr = nullptr;
    // For paper data record
    float alpha_del;  // deleted criterion
    float alpha_bal;  // balance criterion
  };

  struct Operation_Logger_Type {
    PointType point;
    BoxPointType boxpoint;
    bool tree_deleted, tree_downsample_deleted;
    operation_set op;
  };  // 用于记录flatten重新建树时等待的操作

  struct PointType_CMP {
    PointType point;
    float dist = 0.0;
    PointType_CMP(PointType p = PointType(), float d = INFINITY) {
      this->point = p;
      this->dist = d;
    };
    bool operator<(const PointType_CMP& a) const {
      if (fabs(dist - a.dist) < 1e-10)
        return point.x < a.point.x;
      else
        return dist < a.dist;
    }
  };  // 定义一个基于距离的比较器，用于构建最大堆（根节点值最大，父节点值大于子节点的值）

  class MANUAL_HEAP {
   public:
    MANUAL_HEAP(int max_capacity = 100) {
      cap = max_capacity;
      heap = new PointType_CMP[max_capacity];
      heap_size =
          0;  // 构造函数初始化最大堆，但是元素个数仍为0，需要PUSH操作将元素压入堆中
    }

    ~MANUAL_HEAP() { delete[] heap; }

    void pop() {
      if (heap_size == 0) return;
      heap[0] = heap[heap_size -
                     1];  // 用堆中最后一个元素覆盖堆顶元素，即移除堆顶元素
      heap_size--;  // 减少堆中元素个数
      MoveDown(0);  // 从堆顶向下恢复最大堆的性质
      return;
    }

    PointType_CMP top() { return heap[0]; }

    void push(PointType_CMP point) {
      if (heap_size >= cap) return;
      heap[heap_size] = point;  // 将元素压入堆尾
      FloatUp(heap_size);  // 上浮操作，使得当前节点满足最大堆的性质
      heap_size++;
      return;
    }

    int size() { return heap_size; }

    void clear() { heap_size = 0; }

   private:
    int heap_size = 0;
    int cap = 0;
    PointType_CMP* heap;
    // 向下遍历，直至当前节点的值大于或等于其子节点的值 最大堆
    void MoveDown(int heap_index) {
      int l = heap_index * 2 +
              1;  // heap_index为元素再数组表示的堆中索引，l为左子节点索引
      PointType_CMP tmp = heap[heap_index];  // 保存当前节点的值
      while (l < heap_size) {  // 表示当前节点的左子节点存在
        if (l + 1 < heap_size && heap[l] < heap[l + 1])
          l++;  // 右子节点也存在，且左子节点小于右子节点
        if (tmp < heap[l]) {  // 如果父节点小于左子节点
          heap[heap_index] = heap[l];
          heap_index = l;
          l = heap_index * 2 + 1;  // 直到当前节点的值大于或等于其子节点的值
        } else
          break;
      }
      heap[heap_index] = tmp;
      return;
    }
    // 上浮操作，使得当前节点满足最大堆的性质
    void FloatUp(int heap_index) {
      int ancestor = (heap_index - 1) / 2;
      PointType_CMP tmp = heap[heap_index];
      while (heap_index > 0) {
        if (heap[ancestor] < tmp) {
          heap[heap_index] = heap[ancestor];
          heap_index = ancestor;
          ancestor = (heap_index - 1) / 2;
        } else
          break;
      }
      heap[heap_index] = tmp;
      return;
    }
  };

 private:
  // Multi-thread Tree Rebuild
  bool termination_flag = false;
  bool rebuild_flag = false;
  pthread_t rebuild_thread;
  pthread_mutex_t termination_flag_mutex_lock, rebuild_ptr_mutex_lock,
      working_flag_mutex, search_flag_mutex;
  pthread_mutex_t rebuild_logger_mutex_lock, points_deleted_rebuild_mutex_lock;
  // queue<Operation_Logger_Type> Rebuild_Logger;
  MANUAL_Q<Operation_Logger_Type> Rebuild_Logger;  // 重建子树时暂存的增量操作
  PointVector Rebuild_PCL_Storage;  // 保存需要重建的subtree中的点
  KD_TREE_NODE** Rebuild_Ptr = nullptr;  // 指针的指针，指向节点的子节点
  int search_mutex_counter = 0;
  static void* multi_thread_ptr(void* arg);
  void multi_thread_rebuild();
  void start_thread();
  void stop_thread();
  void run_operation(
      KD_TREE_NODE** root,
      Operation_Logger_Type operation);  // 运行暂存的相关增量操作
  // KD Tree Functions and augmented variables
  int Treesize_tmp = 0,
      Validnum_tmp =
          0;  // 临时保存的treesize和validnum，防止多线程时访问出现冲突
  float alpha_bal_tmp = 0.5, alpha_del_tmp = 0.0;
  float delete_criterion_param = 0.5f;
  float balance_criterion_param = 0.7f;
  float downsample_size = 0.2f;
  bool Delete_Storage_Disabled = false;
  KD_TREE_NODE* STATIC_ROOT_NODE = nullptr;  // 虚拟头节点
  PointVector Points_deleted;
  PointVector Downsample_Storage;
  PointVector Multithread_Points_deleted;
  void InitTreeNode(KD_TREE_NODE* root);
  void BuildTree(KD_TREE_NODE** root, int l, int r,
                 PointVector& Storage);  // 构建kdtree
  void Rebuild(KD_TREE_NODE** root);     // 主线程（单线程）重建子树
  int Delete_by_range(KD_TREE_NODE** root, BoxPointType boxpoint,
                      bool allow_rebuild, bool is_downsample);
  void Delete_by_point(KD_TREE_NODE** root, PointType point,
                       bool allow_rebuild);
  void Add_by_point(KD_TREE_NODE** root, PointType point, bool allow_rebuild,
                    int father_axis);
  void Add_by_range(KD_TREE_NODE** root, BoxPointType boxpoint,
                    bool allow_rebuild);
  void Search(KD_TREE_NODE* root, int k_nearest, PointType point,
              MANUAL_HEAP& q,
              double max_dist);  // priority_queue<PointType_CMP>
  void Search_by_range(KD_TREE_NODE* root, BoxPointType boxpoint,
                       PointVector& Storage);
  void Search_by_radius(KD_TREE_NODE* root, PointType point, float radius,
                        PointVector& Storage);
  bool Criterion_Check(KD_TREE_NODE* root);
  void Push_Down(KD_TREE_NODE* root);
  void Update(KD_TREE_NODE* root);
  void delete_tree_nodes(KD_TREE_NODE** root);
  bool same_point(PointType a, PointType b);
  float calc_dist(PointType a, PointType b);
  float calc_box_dist(KD_TREE_NODE* node, PointType point);
  static bool point_cmp_x(PointType a, PointType b);
  static bool point_cmp_y(PointType a, PointType b);
  static bool point_cmp_z(PointType a, PointType b);

 public:
  KD_TREE(float delete_param = 0.5, float balance_param = 0.6,
          float box_length = 0.2);  // 初始化kdtree，将启动multi_rebuild线程
  ~KD_TREE();
  void Set_delete_criterion_param(float delete_param);    // 设置alpha_del
  void Set_balance_criterion_param(float balance_param);  // 设置alpha_bal
  void set_downsample_param(float box_length);  // 设置下采样的分辨率
  void InitializeKDTree(float delete_param = 0.5, float balance_param = 0.7,
                        float box_length = 0.2);  // 执行上面三个函数
  int size();                                     // 返回kdtree的尺寸
  int validnum();  // 返回tree中有效点个数
  void root_alpha(float& alpha_bal,
                  float& alpha_del);    // 获取alpha_bal和alpha_del
  void Build(PointVector point_cloud);  // 构建kdtree，初始建树
  void Nearest_Search(
      PointType point, int k_nearest, PointVector& Nearest_Points,
      vector<float>& Point_Distance,
      double max_dist = INFINITY);  // 最近邻搜索，内部时行search()接口；
                                    // 最大搜索距离默认无穷大
  void Box_Search(const BoxPointType& Box_of_Point,
                  PointVector& Storage);  // 基于BOX搜索，给定BOX的最近最远点；
                                          // 内部执行search_by_range()接口
  void Radius_Search(
      PointType point, const float radius,
      PointVector& Storage);  // 即与半径搜索节点; 内部执行Search_by_radius
  int Add_Points(PointVector& PointToAdd,
                 bool downsample_on);  // 增量式增加点云，需要指定是否进行下采样
  void Add_Point_Boxes(vector<BoxPointType>& BoxPoints);
  void Delete_Points(PointVector& PointToDel);  // 调用delete_by_point接口
  int Delete_Point_Boxes(
      vector<BoxPointType>& BoxPoints);  // 调用delete_by_range
  void flatten(KD_TREE_NODE* root, PointVector& Storage,
               delete_point_storage_set
                   storage_type);  // 遍历tree，将有效节点存入Storage
  void acquire_removed_points(PointVector& removed_points);
  BoxPointType tree_range();  // 统计tree的range
  PointVector PCL_Storage;
  KD_TREE_NODE* Root_Node = nullptr;  // kdtree的根节点
  int max_queue_size = 0;             // 存储Rebuild_Logger队列大小
};
