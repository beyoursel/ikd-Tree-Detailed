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

enum DeletePointStorageSet {
  NOT_RECORD,
  DELETE_POINTS_REC,
};

struct BoxPointType {
  float vertex_min[3];
  float vertex_max[3];
};

template <typename T>
class MANUAL_Q {
 private:
  int head = 0, tail = 0, counter = 0;
  T q[Q_LEN];
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
  using PointVector =
      std::vector<PointType, Eigen::aligned_allocator<PointType>>;
  using Ptr = std::shared_ptr<KD_TREE<PointType>>;
  struct KD_TREE_NODE {
    PointType point;
    uint8_t division_axis;
    int tree_size = 1;
    int invalid_point_num = 0;
    int down_del_num = 0;
    bool point_deleted = false;
    bool tree_deleted = false;
    bool point_downsample_deleted = false;
    bool tree_downsample_deleted = false;
    bool need_push_down_to_left = false;
    bool need_push_down_to_right = false;
    bool working_flag = false;
    float radius_sq;
    pthread_mutex_t push_down_mutex_lock;
    float node_range_x[2], node_range_y[2], node_range_z[2];
    KD_TREE_NODE* left_son_ptr = nullptr;
    KD_TREE_NODE* right_son_ptr = nullptr;
    KD_TREE_NODE* father_ptr = nullptr;
    // For paper data record
    float alpha_del;
    float alpha_bal;
  };

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
  };

  class MANUAL_HEAP {
   public:
    MANUAL_HEAP(int max_capacity = 100) {
      cap = max_capacity;
      heap = new PointType_CMP[max_capacity];
      heap_size = 0;
    }

    ~MANUAL_HEAP() { delete[] heap; }

    void pop() {
      if (heap_size == 0) return;
      heap[0] = heap[heap_size - 1];
      heap_size--;
      MoveDown(0);
      return;
    }

    PointType_CMP top() { return heap[0]; }

    void push(PointType_CMP point) {
      if (heap_size >= cap) return;
      heap[heap_size] = point;
      FloatUp(heap_size);
      heap_size++;
      return;
    }

    int size() { return heap_size; }

    void clear() { heap_size = 0; }

   private:
    int heap_size = 0;
    int cap = 0;
    PointType_CMP* heap;
    void MoveDown(int heap_index) {
      int l = heap_index * 2 + 1;
      PointType_CMP tmp = heap[heap_index];
      while (l < heap_size) {
        if (l + 1 < heap_size && heap[l] < heap[l + 1]) l++;
        if (tmp < heap[l]) {
          heap[heap_index] = heap[l];
          heap_index = l;
          l = heap_index * 2 + 1;
        } else
          break;
      }
      heap[heap_index] = tmp;
      return;
    }

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
  // KD Tree Functions and augmented variables
  float delete_criterion_param_;
  float balance_criterion_param_;
  float downsample_size_;
  KD_TREE_NODE* static_root_node_;
  PointVector points_deleted_;
  PointVector downsample_storage_;
  void InitTreeNode(KD_TREE_NODE* root);
  void BuildTree(KD_TREE_NODE** root, int l, int r, PointVector& Storage);
  void Rebuild(KD_TREE_NODE** root);
  int DeleteByRange(KD_TREE_NODE** root, BoxPointType boxpoint,
                    bool allow_rebuild, bool is_downsample);
  void DeleteByPoint(KD_TREE_NODE** root, PointType point, bool allow_rebuild);
  void AddByPoint(KD_TREE_NODE** root, PointType point, bool allow_rebuild,
                  int father_axis);
  void Search(KD_TREE_NODE* root, int k_nearest, PointType point,
              MANUAL_HEAP& q,
              double max_dist);  // priority_queue<PointType_CMP>
  void SearchByRange(KD_TREE_NODE* root, BoxPointType boxpoint,
                     PointVector& Storage);
  void SearchByRadius(KD_TREE_NODE* root, PointType point, float radius,
                      PointVector& Storage);
  bool CriterionCheck(KD_TREE_NODE* root);
  void PushDown(KD_TREE_NODE* root);
  void Update(KD_TREE_NODE* root);
  void DeleteTreeNodes(KD_TREE_NODE** root);
  bool SamePoint(PointType a, PointType b);
  float CalcDist(PointType a, PointType b);
  float CalcBoxDist(KD_TREE_NODE* node, PointType point);
  static bool PointCmpX(PointType a, PointType b);
  static bool PointCmpY(PointType a, PointType b);
  static bool PointCmpZ(PointType a, PointType b);

 public:
  KD_TREE(float delete_param = 0.5, float balance_param = 0.6,
          float box_length = 0.2);
  ~KD_TREE();
  void SetDeleteCriterionParam(float delete_param);
  void SetBalanceCriterionParam(float balance_param);
  void SetDownsampleParam(float box_length);
  void InitializeKDTree(float delete_param = 0.5, float balance_param = 0.7,
                        float box_length = 0.2);
  int Size();
  int ValidNum();
  void RootAlpha(float& alpha_bal, float& alpha_del);
  void Build(PointVector point_cloud);
  void NearestSearch(PointType point, int k_nearest,
                     PointVector& Nearest_Points,
                     std::vector<float>& Point_Distance,
                     double max_dist = INFINITY);
  void BoxSearch(const BoxPointType& Box_of_Point, PointVector& Storage);
  void RadiusSearch(PointType point, const float radius, PointVector& Storage);
  int AddPoints(PointVector& PointToAdd, bool downsample_on);
  void DeletePoints(PointVector& PointToDel);
  void Flatten(KD_TREE_NODE* root, PointVector& Storage,
               DeletePointStorageSet storage_type);
  void AcquireRemovedPoints(PointVector& removed_points);
  void Clear();
  BoxPointType TreeRange();
  PointVector pcl_storage_;
  KD_TREE_NODE* root_node_;
};
