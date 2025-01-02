#include "ikd_Tree.h"

/*
Description: ikd-Tree: an incremental k-d tree for robotic applications
Author: Yixi Cai
email: yixicai@connect.hku.hk
*/

template <typename PointType>
KD_TREE<PointType>::KD_TREE(float delete_param, float balance_param,
                            float box_length) {
  delete_criterion_param = delete_param;
  balance_criterion_param = balance_param;
  downsample_size = box_length;
  Rebuild_Logger.clear();
  termination_flag = false;
  start_thread(); // 启动multi-rebuild线程，实际上是额外创建了一个线程
}

template <typename PointType> KD_TREE<PointType>::~KD_TREE() {
  stop_thread();
  Delete_Storage_Disabled = true;
  delete_tree_nodes(&Root_Node);
  PointVector().swap(PCL_Storage);
  Rebuild_Logger.clear();
  // malloc_trim(0); // 释放内存
}

template <typename PointType>
void KD_TREE<PointType>::Set_delete_criterion_param(float delete_param) {
  delete_criterion_param = delete_param;
}

template <typename PointType>
void KD_TREE<PointType>::Set_balance_criterion_param(float balance_param) {
  balance_criterion_param = balance_param;
}

template <typename PointType>
void KD_TREE<PointType>::set_downsample_param(float downsample_param) {
  downsample_size = downsample_param;
}

template <typename PointType>
void KD_TREE<PointType>::InitializeKDTree(float delete_param,
                                          float balance_param,
                                          float box_length) {
  Set_delete_criterion_param(delete_param);
  Set_balance_criterion_param(balance_param);
  set_downsample_param(box_length);
}

template <typename PointType>
void KD_TREE<PointType>::InitTreeNode(
    KD_TREE_NODE *root) { // 初始化节点的相关属性
  root->point.x = 0.0f;
  root->point.y = 0.0f;
  root->point.z = 0.0f;
  root->node_range_x[0] = 0.0f;
  root->node_range_x[1] = 0.0f;
  root->node_range_y[0] = 0.0f;
  root->node_range_y[1] = 0.0f;
  root->node_range_z[0] = 0.0f;
  root->node_range_z[1] = 0.0f;
  root->division_axis = 0;
  root->father_ptr = nullptr;
  root->left_son_ptr = nullptr;
  root->right_son_ptr = nullptr;
  root->TreeSize = 0;
  root->invalid_point_num = 0;
  root->down_del_num = 0;
  root->point_deleted = false;
  root->tree_deleted = false;
  root->need_push_down_to_left = false;
  root->need_push_down_to_right = false;
  root->point_downsample_deleted = false;
  root->working_flag = false;
  pthread_mutex_init(&(root->push_down_mutex_lock), NULL);
}

template <typename PointType> int KD_TREE<PointType>::size() {
  int s = 0;
  if (Rebuild_Ptr == nullptr || *Rebuild_Ptr != Root_Node) {
    if (Root_Node != nullptr) {
      return Root_Node->TreeSize;
    } else {
      return 0;
    }
  } else {
    if (!pthread_mutex_trylock(&working_flag_mutex)) {
      s = Root_Node->TreeSize;
      pthread_mutex_unlock(&working_flag_mutex);
      return s;
    } else {
      return Treesize_tmp;
    }
  }
}

template <typename PointType> BoxPointType KD_TREE<PointType>::tree_range() {
  BoxPointType range;
  if (Rebuild_Ptr == nullptr || *Rebuild_Ptr != Root_Node) {
    if (Root_Node != nullptr) {
      range.vertex_min[0] = Root_Node->node_range_x[0];
      range.vertex_min[1] = Root_Node->node_range_y[0];
      range.vertex_min[2] = Root_Node->node_range_z[0];
      range.vertex_max[0] = Root_Node->node_range_x[1];
      range.vertex_max[1] = Root_Node->node_range_y[1];
      range.vertex_max[2] = Root_Node->node_range_z[1];
    } else {
      memset(&range, 0, sizeof(range));
    }
  } else {
    if (!pthread_mutex_trylock(&working_flag_mutex)) {
      range.vertex_min[0] = Root_Node->node_range_x[0];
      range.vertex_min[1] = Root_Node->node_range_y[0];
      range.vertex_min[2] = Root_Node->node_range_z[0];
      range.vertex_max[0] = Root_Node->node_range_x[1];
      range.vertex_max[1] = Root_Node->node_range_y[1];
      range.vertex_max[2] = Root_Node->node_range_z[1];
      pthread_mutex_unlock(&working_flag_mutex);
    } else {
      memset(&range, 0, sizeof(range));
    }
  }
  return range;
}

template <typename PointType> int KD_TREE<PointType>::validnum() {
  int s = 0;
  if (Rebuild_Ptr == nullptr || *Rebuild_Ptr != Root_Node) {
    if (Root_Node != nullptr)
      return (Root_Node->TreeSize - Root_Node->invalid_point_num);
    else
      return 0;
  } else {
    if (!pthread_mutex_trylock(&working_flag_mutex)) {
      s = Root_Node->TreeSize - Root_Node->invalid_point_num;
      pthread_mutex_unlock(&working_flag_mutex);
      return s;
    } else {
      return -1;
    }
  }
}

template <typename PointType>
void KD_TREE<PointType>::root_alpha(float &alpha_bal, float &alpha_del) {
  if (Rebuild_Ptr == nullptr ||
      *Rebuild_Ptr != Root_Node) { // 重建节点不存在或者不为Root_Node
    alpha_bal = Root_Node->alpha_bal;
    alpha_del = Root_Node->alpha_del;
    return;
  } else {
    if (!pthread_mutex_trylock(
            &working_flag_mutex)) { // 非阻塞式加锁，加锁成功返回0
      alpha_bal = Root_Node->alpha_bal;
      alpha_del = Root_Node->alpha_del;
      pthread_mutex_unlock(&working_flag_mutex);
      return;
    } else {
      alpha_bal = alpha_bal_tmp;
      alpha_del = alpha_del_tmp;
      return;
    }
  }
}

template <typename PointType>
void KD_TREE<PointType>::start_thread() { // 初始化互斥锁，创建线程
  pthread_mutex_init(&termination_flag_mutex_lock, NULL);
  pthread_mutex_init(&rebuild_ptr_mutex_lock, NULL);
  pthread_mutex_init(&rebuild_logger_mutex_lock, NULL);
  pthread_mutex_init(&points_deleted_rebuild_mutex_lock, NULL);
  pthread_mutex_init(&working_flag_mutex, NULL);
  pthread_mutex_init(&search_flag_mutex, NULL);
  pthread_create(
      &rebuild_thread, NULL, multi_thread_ptr,
      (void *)this); // multi_thread_ptr为线程入口函数，即线程需要执行的代码逻辑
  printf("Multi thread started \n");
}

template <typename PointType> void KD_TREE<PointType>::stop_thread() {
  pthread_mutex_lock(&termination_flag_mutex_lock);
  termination_flag = true;
  pthread_mutex_unlock(&termination_flag_mutex_lock);
  if (rebuild_thread)
    pthread_join(rebuild_thread, NULL);
  pthread_mutex_destroy(&termination_flag_mutex_lock);
  pthread_mutex_destroy(&rebuild_logger_mutex_lock);
  pthread_mutex_destroy(&rebuild_ptr_mutex_lock);
  pthread_mutex_destroy(&points_deleted_rebuild_mutex_lock);
  pthread_mutex_destroy(&working_flag_mutex);
  pthread_mutex_destroy(&search_flag_mutex);
}

template <typename PointType>
void *KD_TREE<PointType>::multi_thread_ptr(void *arg) {
  KD_TREE *handle = (KD_TREE *)arg;
  handle->multi_thread_rebuild();
  return nullptr;
}

template <typename PointType> void KD_TREE<PointType>::multi_thread_rebuild() {
  /**
   * \brief
   * 当rebuilding比较大的subtree时，耗时严重，因此需要采用多线程重建，当subtree的节点数量少于Nmax时，使用单线程重建，即主线程，否则启动多线程（对应论文中的second
   * thread） 0. 多线程重建时会lock除去query的其他增量操作（points insert,
   * re-insert, delete）
   * 1.
   * 首先copy子树中所有有效节点到Vector中(即flatten)，在rebuilding过程中保持原来的子树不变。
   * 2.
   * flatten过后，subtree被解锁，继续执行增量操作（按论文中所说时在原来的subtree进行，实际上没有，仅在平衡的subtree上执行）,这些未执行的操作暂存在Rebuild_Logger。
   * 3. 一但平衡后的subtree构建完毕，将Rebuild_Logger中的操作依次执行。
   */
  bool terminated = false;
  KD_TREE_NODE *father_ptr, **new_node_ptr;
  pthread_mutex_lock(&termination_flag_mutex_lock);
  terminated = termination_flag; // 根据termination_flag判断是否rebuild
  pthread_mutex_unlock(&termination_flag_mutex_lock);
  while (!terminated) {
    pthread_mutex_lock(&rebuild_ptr_mutex_lock);
    pthread_mutex_lock(&working_flag_mutex);
    if (Rebuild_Ptr !=
        nullptr) { // 当Rebuild_Ptr不为空时，表示需要重建子树，Rebuild_Ptr是在主线程中赋值的
      /* Traverse and copy */
      if (!Rebuild_Logger.empty()) {
        printf("\n\n\n\n\n\n\n\n\n\n\n ERROR!!! \n\n\n\n\n\n\n\n\n");
      }
      rebuild_flag = true;
      if (*Rebuild_Ptr ==
          Root_Node) { // 当Rebuild_Ptr指向root
                       // node时，重建树的相关属性直接使用root_node的
        Treesize_tmp = Root_Node->TreeSize;
        Validnum_tmp = Root_Node->TreeSize - Root_Node->invalid_point_num;
        alpha_bal_tmp = Root_Node->alpha_bal;
        alpha_del_tmp = Root_Node->alpha_del;
      }
      KD_TREE_NODE *old_root_node = (*Rebuild_Ptr); // 存储原始根节点
      father_ptr = (*Rebuild_Ptr)->father_ptr; // 存储原始根节点的父节点
      PointVector().swap(
          Rebuild_PCL_Storage); // 释放Rebuild_PCL_Storage中的内存
      // Lock Search
      pthread_mutex_lock(&search_flag_mutex);
      while (search_mutex_counter != 0) { // 当不为0时，表示别的线程正在占用资源
        pthread_mutex_unlock(&search_flag_mutex); //临时解锁
        usleep(1);                                // 休眠1 us
        pthread_mutex_lock(
            &search_flag_mutex); // 加锁，保证读取search_mutex_counter时安全
      }
      search_mutex_counter = -1;
      pthread_mutex_unlock(&search_flag_mutex);
      // Lock deleted points cache
      pthread_mutex_lock(
          &points_deleted_rebuild_mutex_lock); // flatten操作需要加锁
      flatten(
          *Rebuild_Ptr, Rebuild_PCL_Storage,
          MULTI_THREAD_REC); // 将Rebuild_Ptr的节点按序存入Rebuild_PCL_Storage
      // Unlock deleted points cache
      pthread_mutex_unlock(&points_deleted_rebuild_mutex_lock);
      // Unlock Search
      pthread_mutex_lock(&search_flag_mutex); // 多线程访问共享资源时需要加锁
      search_mutex_counter = 0;
      pthread_mutex_unlock(&search_flag_mutex);
      pthread_mutex_unlock(&working_flag_mutex);
      /* Rebuild and update missed operations*/
      Operation_Logger_Type Operation;
      KD_TREE_NODE *new_root_node = nullptr;     // 子树新的根节点
      if (int(Rebuild_PCL_Storage.size()) > 0) { // 若满足重新建树条件
        BuildTree(&new_root_node, 0, Rebuild_PCL_Storage.size() - 1,
                  Rebuild_PCL_Storage); // 构建树
        // Rebuild has been done. Updates the blocked operations into the new
        // tree
        pthread_mutex_lock(&working_flag_mutex);
        pthread_mutex_lock(&rebuild_logger_mutex_lock);
        int tmp_counter = 0;
        while (!Rebuild_Logger.empty()) { // 遍历Rebuild_Logger中存储的操作
          Operation = Rebuild_Logger.front();
          max_queue_size = max(max_queue_size, Rebuild_Logger.size());
          Rebuild_Logger.pop();
          pthread_mutex_unlock(&rebuild_logger_mutex_lock);
          pthread_mutex_unlock(&working_flag_mutex);
          run_operation(&new_root_node,
                        Operation); // 在重建后的子树上执行队列中的操作
          tmp_counter++;
          if (tmp_counter % 10 == 0)
            usleep(1); // 每10次休眠1us
          pthread_mutex_lock(&working_flag_mutex);
          pthread_mutex_lock(&rebuild_logger_mutex_lock);
        }
        pthread_mutex_unlock(&rebuild_logger_mutex_lock);
      }
      /* Replace to original tree*/
      // pthread_mutex_lock(&working_flag_mutex);
      pthread_mutex_lock(&search_flag_mutex);
      while (search_mutex_counter != 0) {
        pthread_mutex_unlock(&search_flag_mutex);
        usleep(1);
        pthread_mutex_lock(&search_flag_mutex);
      }
      search_mutex_counter = -1;
      pthread_mutex_unlock(&search_flag_mutex);
      if (father_ptr->left_son_ptr ==
          *Rebuild_Ptr) { // 将重建的subtree替代原来的subtree
        father_ptr->left_son_ptr = new_root_node;
      } else if (father_ptr->right_son_ptr == *Rebuild_Ptr) {
        father_ptr->right_son_ptr = new_root_node;
      } else {
        throw "Error: Father ptr incompatible with current node\n";
      }
      if (new_root_node != nullptr)
        new_root_node->father_ptr = father_ptr;
      (*Rebuild_Ptr) = new_root_node;
      int valid_old =
          old_root_node->TreeSize - old_root_node->invalid_point_num;
      int valid_new =
          new_root_node->TreeSize - new_root_node->invalid_point_num;
      if (father_ptr == STATIC_ROOT_NODE)
        Root_Node = STATIC_ROOT_NODE->left_son_ptr; // ？？？
      KD_TREE_NODE *update_root = *Rebuild_Ptr;
      while (
          update_root != nullptr &&
          update_root !=
              Root_Node) { // 向上遍历更新各个节点的信息,range、tree_size、alpha_bal、alpha_del等
        update_root = update_root->father_ptr;
        if (update_root->working_flag)
          break;
        if (update_root == update_root->father_ptr->left_son_ptr &&
            update_root->father_ptr->need_push_down_to_left)
          break;
        if (update_root == update_root->father_ptr->right_son_ptr &&
            update_root->father_ptr->need_push_down_to_right)
          break;
        Update(update_root); // 更新节点属性
      }
      pthread_mutex_lock(&search_flag_mutex);
      search_mutex_counter = 0;
      pthread_mutex_unlock(&search_flag_mutex);
      Rebuild_Ptr = nullptr;
      pthread_mutex_unlock(&working_flag_mutex);
      rebuild_flag = false;
      /* Delete discarded tree nodes */
      delete_tree_nodes(
          &old_root_node); // 完成子树重建之后，删除deleted标记为true的节点
    } else {
      pthread_mutex_unlock(&working_flag_mutex);
    }
    pthread_mutex_unlock(&rebuild_ptr_mutex_lock);
    pthread_mutex_lock(&termination_flag_mutex_lock);
    terminated = termination_flag;
    pthread_mutex_unlock(&termination_flag_mutex_lock);
    usleep(100);
  }
  printf("Rebuild thread terminated normally\n");
}

template <typename PointType>
void KD_TREE<PointType>::run_operation(KD_TREE_NODE **root,
                                       Operation_Logger_Type operation) {
  switch (operation.op) {
  case ADD_POINT:
    Add_by_point(root, operation.point, false, (*root)->division_axis);
    break;
  case ADD_BOX:
    Add_by_range(root, operation.boxpoint, false);
    break;
  case DELETE_POINT:
    Delete_by_point(root, operation.point, false);
    break;
  case DELETE_BOX:
    Delete_by_range(root, operation.boxpoint, false, false);
    break;
  case DOWNSAMPLE_DELETE:
    Delete_by_range(root, operation.boxpoint, false, true);
    break;
  case PUSH_DOWN:
    (*root)->tree_downsample_deleted |= operation.tree_downsample_deleted;
    (*root)->point_downsample_deleted |= operation.tree_downsample_deleted;
    (*root)->tree_deleted =
        operation.tree_deleted || (*root)->tree_downsample_deleted;
    (*root)->point_deleted =
        (*root)->tree_deleted || (*root)->point_downsample_deleted;
    if (operation.tree_downsample_deleted)
      (*root)->down_del_num = (*root)->TreeSize;
    if (operation.tree_deleted)
      (*root)->invalid_point_num = (*root)->TreeSize;
    else
      (*root)->invalid_point_num = (*root)->down_del_num;
    (*root)->need_push_down_to_left = true;
    (*root)->need_push_down_to_right = true;
    break;
  default:
    break;
  }
}

template <typename PointType>
void KD_TREE<PointType>::Build(PointVector point_cloud) {
  /**
   * 1. 清空之前存在的kdtree
   * 2. 建立虚拟头节点并初始化
   * 3. 根据输入点云建立kdtree
   * 4. 更新虚拟节点的相关属性
   * 5. 虚拟节点左孩子指向根节点
   */
  if (Root_Node != nullptr) { // 初始创建树时，需要清空之前的树
    delete_tree_nodes(&Root_Node);
  } // 若kdtree已经存在，则从Root_Node开始清空
  if (point_cloud.size() == 0)
    return;
  STATIC_ROOT_NODE = new KD_TREE_NODE; // 创建虚拟头节点
  InitTreeNode(STATIC_ROOT_NODE);      // 初始化虚拟头节点
  BuildTree(&STATIC_ROOT_NODE->left_son_ptr, 0, point_cloud.size() - 1,
            point_cloud);         // build-tree
  Update(STATIC_ROOT_NODE);       // 更新节点属性
  STATIC_ROOT_NODE->TreeSize = 0; // 虚拟头节点的treesize设置为0
  Root_Node =
      STATIC_ROOT_NODE->left_son_ptr; // kdtree的根节点设置为虚拟头节点的左孩子
}

template <typename PointType>
void KD_TREE<PointType>::Nearest_Search(PointType point, int k_nearest,
                                        PointVector &Nearest_Points,
                                        vector<float> &Point_Distance,
                                        double max_dist) {
  MANUAL_HEAP q(2 * k_nearest); // 维护一个优先级队列，最大堆
  q.clear();
  vector<float>().swap(Point_Distance);
  if (Rebuild_Ptr == nullptr || *Rebuild_Ptr != Root_Node) {
    Search(Root_Node, k_nearest, point, q, max_dist);
  } else {
    pthread_mutex_lock(&search_flag_mutex);
    while (search_mutex_counter == -1) {
      pthread_mutex_unlock(&search_flag_mutex);
      usleep(1);
      pthread_mutex_lock(&search_flag_mutex);
    }
    search_mutex_counter += 1;
    pthread_mutex_unlock(&search_flag_mutex);
    Search(Root_Node, k_nearest, point, q, max_dist);
    pthread_mutex_lock(&search_flag_mutex);
    search_mutex_counter -= 1;
    pthread_mutex_unlock(&search_flag_mutex);
  }
  int k_found = min(k_nearest, int(q.size()));
  PointVector().swap(Nearest_Points);
  vector<float>().swap(Point_Distance);
  for (int i = 0; i < k_found; i++) {
    Nearest_Points.insert(Nearest_Points.begin(), q.top().point);
    Point_Distance.insert(Point_Distance.begin(), q.top().dist);
    q.pop();
  }
  return;
}

template <typename PointType>
void KD_TREE<PointType>::Box_Search(const BoxPointType &Box_of_Point,
                                    PointVector &Storage) {
  Storage.clear();
  Search_by_range(Root_Node, Box_of_Point, Storage);
}

template <typename PointType>
void KD_TREE<PointType>::Radius_Search(PointType point, const float radius,
                                       PointVector &Storage) {
  Storage.clear();
  Search_by_radius(Root_Node, point, radius, Storage);
}

template <typename PointType>
int KD_TREE<PointType>::Add_Points(PointVector &PointToAdd,
                                   bool downsample_on) {
  /**
   * \brief 向kdtree中增加点云，遍历输入点云，调用Add_By_Point
   * 1. 是否进行下采样
   * 2.
   * 若进行下采样，根据降采样分辨率，找到插入点对应的bbox，然后计算bbox的中点，搜索bbox内的点云，计算距离中点最近的点，作为降采样后的点。
   * 3. 删除降采样后的无效点
   */
  int NewPointSize = PointToAdd.size();
  int tree_size = size(); // 获取当前kdtree的节点个数
  BoxPointType Box_of_Point;
  PointType downsample_result, mid_point;
  bool downsample_switch = downsample_on && DOWNSAMPLE_SWITCH;
  float min_dist, tmp_dist;
  int tmp_counter = 0;
  for (int i = 0; i < PointToAdd.size(); i++) { // 遍历新增的点云
    if (downsample_switch) {                    // 若需要下采样
      // 计算当前点落入的box
      Box_of_Point.vertex_min[0] =
          floor(PointToAdd[i].x / downsample_size) * downsample_size;
      Box_of_Point.vertex_max[0] = Box_of_Point.vertex_min[0] + downsample_size;
      Box_of_Point.vertex_min[1] =
          floor(PointToAdd[i].y / downsample_size) * downsample_size;
      Box_of_Point.vertex_max[1] = Box_of_Point.vertex_min[1] + downsample_size;
      Box_of_Point.vertex_min[2] =
          floor(PointToAdd[i].z / downsample_size) * downsample_size;
      Box_of_Point.vertex_max[2] = Box_of_Point.vertex_min[2] + downsample_size;
      // 计算box中心点
      mid_point.x =
          Box_of_Point.vertex_min[0] +
          (Box_of_Point.vertex_max[0] - Box_of_Point.vertex_min[0]) / 2.0;
      mid_point.y =
          Box_of_Point.vertex_min[1] +
          (Box_of_Point.vertex_max[1] - Box_of_Point.vertex_min[1]) / 2.0;
      mid_point.z =
          Box_of_Point.vertex_min[2] +
          (Box_of_Point.vertex_max[2] - Box_of_Point.vertex_min[2]) / 2.0;
      PointVector().swap(Downsample_Storage);
      // 搜索box内的节点坐标
      Search_by_range(Root_Node, Box_of_Point, Downsample_Storage);
      min_dist =
          calc_dist(PointToAdd[i], mid_point); // 计算当前增加点与box中点的距离
      downsample_result = PointToAdd[i];
      for (
          int index = 0; index < Downsample_Storage.size();
          index++) { // 遍历box中的点，选取距离mid_point最近的点作为下采样后的点
        tmp_dist = calc_dist(Downsample_Storage[index], mid_point);
        if (tmp_dist < min_dist) {
          min_dist = tmp_dist;
          downsample_result = Downsample_Storage[index];
        }
      }
      if (Rebuild_Ptr == nullptr || *Rebuild_Ptr != Root_Node) {
        if (Downsample_Storage.size() > 1 ||
            same_point(PointToAdd[i], downsample_result)) {
          if (Downsample_Storage.size() > 0)
            Delete_by_range(&Root_Node, Box_of_Point, true,
                            true); // 删除下采样后不需要的点
          Add_by_point(
              &Root_Node, downsample_result, true,
              Root_Node
                  ->division_axis); // 若没有启用多线程重建，这里allow_rebuild为true
          tmp_counter++;
        }
      } else {
        if (Downsample_Storage.size() > 1 ||
            same_point(PointToAdd[i], downsample_result)) {
          Operation_Logger_Type operation_delete, operation;
          operation_delete.boxpoint = Box_of_Point;
          operation_delete.op = DOWNSAMPLE_DELETE;
          operation.point = downsample_result;
          operation.op = ADD_POINT;
          pthread_mutex_lock(&working_flag_mutex);
          if (Downsample_Storage.size() > 0)
            Delete_by_range(&Root_Node, Box_of_Point, false, true);
          Add_by_point(&Root_Node, downsample_result, false,
                       Root_Node->division_axis);
          tmp_counter++;
          if (rebuild_flag) {
            pthread_mutex_lock(&rebuild_logger_mutex_lock);
            if (Downsample_Storage.size() > 0)
              Rebuild_Logger.push(operation_delete);
            Rebuild_Logger.push(operation);
            pthread_mutex_unlock(&rebuild_logger_mutex_lock);
          }
          pthread_mutex_unlock(&working_flag_mutex);
        };
      }
    } else { // 不需要下采样
      if (Rebuild_Ptr == nullptr || *Rebuild_Ptr != Root_Node) {
        Add_by_point(&Root_Node, PointToAdd[i], true, Root_Node->division_axis);
      } else {
        Operation_Logger_Type operation;
        operation.point = PointToAdd[i];
        operation.op = ADD_POINT;
        pthread_mutex_lock(&working_flag_mutex);
        Add_by_point(&Root_Node, PointToAdd[i], false,
                     Root_Node->division_axis);
        if (rebuild_flag) {
          pthread_mutex_lock(&rebuild_logger_mutex_lock);
          Rebuild_Logger.push(operation);
          pthread_mutex_unlock(&rebuild_logger_mutex_lock);
        }
        pthread_mutex_unlock(&working_flag_mutex);
      }
    }
  }
  return tmp_counter;
}

template <typename PointType>
void KD_TREE<PointType>::Add_Point_Boxes(vector<BoxPointType> &BoxPoints) {
  for (int i = 0; i < BoxPoints.size(); i++) {
    if (Rebuild_Ptr == nullptr || *Rebuild_Ptr != Root_Node) {
      Add_by_range(&Root_Node, BoxPoints[i], true);
    } else {
      Operation_Logger_Type operation;
      operation.boxpoint = BoxPoints[i];
      operation.op = ADD_BOX;
      pthread_mutex_lock(&working_flag_mutex);
      Add_by_range(&Root_Node, BoxPoints[i], false);
      if (rebuild_flag) {
        pthread_mutex_lock(&rebuild_logger_mutex_lock);
        Rebuild_Logger.push(operation);
        pthread_mutex_unlock(&rebuild_logger_mutex_lock);
      }
      pthread_mutex_unlock(&working_flag_mutex);
    }
  }
  return;
}

template <typename PointType>
void KD_TREE<PointType>::Delete_Points(PointVector &PointToDel) {
  for (int i = 0; i < PointToDel.size(); i++) {
    if (Rebuild_Ptr == nullptr || *Rebuild_Ptr != Root_Node) {
      Delete_by_point(&Root_Node, PointToDel[i], true);
    } else {
      Operation_Logger_Type operation;
      operation.point = PointToDel[i];
      operation.op = DELETE_POINT;
      pthread_mutex_lock(&working_flag_mutex);
      Delete_by_point(&Root_Node, PointToDel[i], false);
      if (rebuild_flag) {
        pthread_mutex_lock(&rebuild_logger_mutex_lock);
        Rebuild_Logger.push(operation);
        pthread_mutex_unlock(&rebuild_logger_mutex_lock);
      }
      pthread_mutex_unlock(&working_flag_mutex);
    }
  }
  return;
}

template <typename PointType>
int KD_TREE<PointType>::Delete_Point_Boxes(vector<BoxPointType> &BoxPoints) {
  int tmp_counter = 0;
  for (int i = 0; i < BoxPoints.size(); i++) {
    if (Rebuild_Ptr == nullptr || *Rebuild_Ptr != Root_Node) {
      tmp_counter += Delete_by_range(&Root_Node, BoxPoints[i], true, false); // allowed_rebuild = true, downsample = false
    } else {
      Operation_Logger_Type operation;
      operation.boxpoint = BoxPoints[i];
      operation.op = DELETE_BOX;
      pthread_mutex_lock(&working_flag_mutex);
      tmp_counter += Delete_by_range(&Root_Node, BoxPoints[i], false, false);
      if (rebuild_flag) {
        pthread_mutex_lock(&rebuild_logger_mutex_lock);
        Rebuild_Logger.push(operation);
        pthread_mutex_unlock(&rebuild_logger_mutex_lock);
      }
      pthread_mutex_unlock(&working_flag_mutex);
    }
  }
  return tmp_counter;
}

template <typename PointType>
void KD_TREE<PointType>::acquire_removed_points(PointVector &removed_points) {
  // 将points_deleted和multithread_points_deleted中被删除的点合并
  pthread_mutex_lock(&points_deleted_rebuild_mutex_lock);
  for (int i = 0; i < Points_deleted.size(); i++) {
    removed_points.push_back(Points_deleted[i]);
  }
  for (int i = 0; i < Multithread_Points_deleted.size(); i++) {
    removed_points.push_back(Multithread_Points_deleted[i]);
  }
  Points_deleted.clear();
  Multithread_Points_deleted.clear();
  pthread_mutex_unlock(&points_deleted_rebuild_mutex_lock);
  return;
}

template <typename PointType>
void KD_TREE<PointType>::BuildTree(KD_TREE_NODE **root, int l, int r,
                                   PointVector &Storage) {
  if (l > r)
    return;                 // l和r分别为指向Storage头、尾索引
  *root = new KD_TREE_NODE; // 分配内存
  InitTreeNode(*root);      // 初始化根节点相关属性
  int mid = (l + r) >> 1; // 位运算，向右移一位，相当于/，比/高效
  int div_axis = 0;
  int i;
  // Find the best division Axis
  float min_value[3] = {INFINITY, INFINITY, INFINITY};
  float max_value[3] = {-INFINITY, -INFINITY, -INFINITY};
  float dim_range[3] = {0, 0, 0};
  for (i = l; i <= r; i++) { // 统计输入点云在x, y, z三个方向的范围
    min_value[0] = min(min_value[0], Storage[i].x);
    min_value[1] = min(min_value[1], Storage[i].y);
    min_value[2] = min(min_value[2], Storage[i].z);
    max_value[0] = max(max_value[0], Storage[i].x);
    max_value[1] = max(max_value[1], Storage[i].y);
    max_value[2] = max(max_value[2], Storage[i].z);
  }
  // Select the longest dimension as division axis
  for (i = 0; i < 3; i++)
    dim_range[i] = max_value[i] - min_value[i];
  for (i = 1; i < 3; i++)
    if (dim_range[i] > dim_range[div_axis])
      div_axis = i; // 选取点云分布最广的维度为划分轴
  // Divide by the division axis and recursively build.

  (*root)->division_axis = div_axis;
  switch (div_axis) {
  case 0:
    nth_element(begin(Storage) + l, begin(Storage) + mid,
                begin(Storage) + r + 1, point_cmp_x);
    break;
  case 1:
    nth_element(begin(Storage) + l, begin(Storage) + mid,
                begin(Storage) + r + 1, point_cmp_y);
    break;
  case 2:
    nth_element(begin(Storage) + l, begin(Storage) + mid,
                begin(Storage) + r + 1, point_cmp_z);
    break;
  default:
    nth_element(begin(Storage) + l, begin(Storage) + mid,
                begin(Storage) + r + 1, point_cmp_x);
    break;
  } // nth_element()函数可以保证Storage[mid]左侧元素比自身小，而右侧比自身大
  (*root)->point = Storage[mid]; // 将中位数作为根节点
  KD_TREE_NODE *left_son = nullptr, *right_son = nullptr;
  BuildTree(&left_son, l, mid - 1, Storage);  // 递归构建左子树
  BuildTree(&right_son, mid + 1, r, Storage); // 递归构建右子树
  (*root)->left_son_ptr = left_son;
  (*root)->right_son_ptr = right_son;
  Update((
      *root)); // 更新subtree的相关属性，如range、alpha_del、alpha_bal、treesize等
  return;
}

template <typename PointType>
void KD_TREE<PointType>::Rebuild(KD_TREE_NODE **root) {
  KD_TREE_NODE *father_ptr;
  if ((*root)->TreeSize >= Multi_Thread_Rebuild_Point_Num) {
    if (!pthread_mutex_trylock(&rebuild_ptr_mutex_lock)) {
      if (Rebuild_Ptr == nullptr ||
          ((*root)->TreeSize >
           (*Rebuild_Ptr)
               ->TreeSize)) { // 若Rebuild_Ptr为空，或者当前节点的treesize大于当前重建节点的treesize，则将当前节点设为重建节点
        Rebuild_Ptr = root;
      }
      pthread_mutex_unlock(&rebuild_ptr_mutex_lock);
    }
  } else {
    father_ptr = (*root)->father_ptr;
    int size_rec = (*root)->TreeSize;
    PCL_Storage.clear();
    flatten(*root, PCL_Storage, DELETE_POINTS_REC);
    delete_tree_nodes(root);
    BuildTree(root, 0, PCL_Storage.size() - 1, PCL_Storage);
    if (*root != nullptr)
      (*root)->father_ptr = father_ptr;
    if (*root == Root_Node)
      STATIC_ROOT_NODE->left_son_ptr = *root; // STATIC_ROOT_NODE类似虚拟头节点
  }
  return;
}

template <typename PointType>
int KD_TREE<PointType>::Delete_by_range(KD_TREE_NODE **root,
                                        BoxPointType boxpoint,
                                        bool allow_rebuild,
                                        bool is_downsample) {
  if ((*root) == nullptr || (*root)->tree_deleted)
    return 0;
  (*root)->working_flag = true;
  Push_Down(*root);
  int tmp_counter = 0; // 记录被删除的节点数量
  if (boxpoint.vertex_max[0] <= (*root)->node_range_x[0] ||
      boxpoint.vertex_min[0] > (*root)->node_range_x[1])
    return 0;
  if (boxpoint.vertex_max[1] <= (*root)->node_range_y[0] ||
      boxpoint.vertex_min[1] > (*root)->node_range_y[1])
    return 0;
  if (boxpoint.vertex_max[2] <= (*root)->node_range_z[0] ||
      boxpoint.vertex_min[2] > (*root)->node_range_z[1])
    return 0;
  if (boxpoint.vertex_min[0] <= (*root)->node_range_x[0] &&
      boxpoint.vertex_max[0] > (*root)->node_range_x[1] &&
      boxpoint.vertex_min[1] <= (*root)->node_range_y[0] &&
      boxpoint.vertex_max[1] > (*root)->node_range_y[1] &&
      boxpoint.vertex_min[2] <= (*root)->node_range_z[0] &&
      boxpoint.vertex_max[2] > (*root)->node_range_z[1]) {
    // subtree含于box； 这种清空便体现了push_down的作用
    (*root)->tree_deleted = true;
    (*root)->point_deleted = true;
    (*root)->need_push_down_to_left = true;
    (*root)->need_push_down_to_right = true;
    tmp_counter = (*root)->TreeSize - (*root)->invalid_point_num;
    (*root)->invalid_point_num = (*root)->TreeSize;
    if (is_downsample) {
      (*root)->tree_downsample_deleted = true;
      (*root)->point_downsample_deleted = true;
      (*root)->down_del_num = (*root)->TreeSize;
    }
    return tmp_counter;
  }
  if (!(*root)->point_deleted && boxpoint.vertex_min[0] <= (*root)->point.x &&
      boxpoint.vertex_max[0] > (*root)->point.x &&
      boxpoint.vertex_min[1] <= (*root)->point.y &&
      boxpoint.vertex_max[1] > (*root)->point.y &&
      boxpoint.vertex_min[2] <= (*root)->point.z &&
      boxpoint.vertex_max[2] > (*root)->point.z) { // 当前节点在box内
    (*root)->point_deleted = true;
    tmp_counter += 1;
    if (is_downsample)
      (*root)->point_downsample_deleted = true;
  } // 删除当前节点坐标，delete_by_range最小单元
  Operation_Logger_Type delete_box_log;
  struct timespec Timeout;
  if (is_downsample)
    delete_box_log.op = DOWNSAMPLE_DELETE;
  else
    delete_box_log.op = DELETE_BOX;
  delete_box_log.boxpoint = boxpoint;
  // 可以发现，delete_by_range时左右子树都要遍历
  if ((Rebuild_Ptr == nullptr) ||
      (*root)->left_son_ptr != *Rebuild_Ptr) { // 左子树
    tmp_counter += Delete_by_range(&((*root)->left_son_ptr), boxpoint,
                                   allow_rebuild, is_downsample);
  } else {
    pthread_mutex_lock(&working_flag_mutex);
    tmp_counter += Delete_by_range(&((*root)->left_son_ptr), boxpoint, false,
                                   is_downsample);
    if (rebuild_flag) {
      pthread_mutex_lock(&rebuild_logger_mutex_lock);
      Rebuild_Logger.push(delete_box_log);
      pthread_mutex_unlock(&rebuild_logger_mutex_lock);
    }
    pthread_mutex_unlock(&working_flag_mutex);
  }
  if ((Rebuild_Ptr == nullptr) ||
      (*root)->right_son_ptr != *Rebuild_Ptr) { // 右子树
    tmp_counter += Delete_by_range(&((*root)->right_son_ptr), boxpoint,
                                   allow_rebuild, is_downsample);
  } else {
    pthread_mutex_lock(&working_flag_mutex);
    tmp_counter += Delete_by_range(&((*root)->right_son_ptr), boxpoint, false,
                                   is_downsample);
    if (rebuild_flag) {
      pthread_mutex_lock(&rebuild_logger_mutex_lock);
      Rebuild_Logger.push(delete_box_log);
      pthread_mutex_unlock(&rebuild_logger_mutex_lock);
    }
    pthread_mutex_unlock(&working_flag_mutex);
  }
  Update(*root);
  if (Rebuild_Ptr != nullptr && *Rebuild_Ptr == *root &&
      (*root)->TreeSize < Multi_Thread_Rebuild_Point_Num)
    Rebuild_Ptr = nullptr;
  bool need_rebuild = allow_rebuild & Criterion_Check((*root));
  if (need_rebuild)
    Rebuild(root); // 判断是否需要重新建树
  if ((*root) != nullptr)
    (*root)->working_flag = false;
  return tmp_counter;
}

template <typename PointType>
void KD_TREE<PointType>::Delete_by_point(KD_TREE_NODE **root, PointType point,
                                         bool allow_rebuild) {
  if ((*root) == nullptr || (*root)->tree_deleted)
    return; // 当root为空或者该节点下的subtree被删除
  (*root)->working_flag = true;
  Push_Down(*root);
  if (same_point((*root)->point, point) && !(*root)->point_deleted) { // 待删除点为当前节点，且当前节点未删除
    (*root)->point_deleted = true;
    (*root)->invalid_point_num += 1;
    if ((*root)->invalid_point_num == (*root)->TreeSize)
      (*root)->tree_deleted = true;
    return;
  }
  Operation_Logger_Type delete_log;
  struct timespec Timeout;
  delete_log.op = DELETE_POINT;
  delete_log.point = point;
  if (((*root)->division_axis == 0 && point.x < (*root)->point.x) ||
      ((*root)->division_axis == 1 && point.y < (*root)->point.y) ||
      ((*root)->division_axis == 2 && point.z < (*root)->point.z)) { // 左子树范围
    if ((Rebuild_Ptr == nullptr) || (*root)->left_son_ptr != *Rebuild_Ptr) {
      Delete_by_point(&(*root)->left_son_ptr, point, allow_rebuild);
    } else {
      pthread_mutex_lock(&working_flag_mutex);
      Delete_by_point(&(*root)->left_son_ptr, point, false);
      if (rebuild_flag) {
        pthread_mutex_lock(&rebuild_logger_mutex_lock);
        Rebuild_Logger.push(delete_log);
        pthread_mutex_unlock(&rebuild_logger_mutex_lock);
      }
      pthread_mutex_unlock(&working_flag_mutex);
    }
  } else { // 右子树范围
    if ((Rebuild_Ptr == nullptr) || (*root)->right_son_ptr != *Rebuild_Ptr) {
      Delete_by_point(&(*root)->right_son_ptr, point, allow_rebuild);
    } else {
      pthread_mutex_lock(&working_flag_mutex);
      Delete_by_point(&(*root)->right_son_ptr, point, false);
      if (rebuild_flag) {
        pthread_mutex_lock(&rebuild_logger_mutex_lock);
        Rebuild_Logger.push(delete_log);
        pthread_mutex_unlock(&rebuild_logger_mutex_lock);
      }
      pthread_mutex_unlock(&working_flag_mutex);
    }
  }
  Update(*root);
  if (Rebuild_Ptr != nullptr && *Rebuild_Ptr == *root &&
      (*root)->TreeSize < Multi_Thread_Rebuild_Point_Num)
    Rebuild_Ptr = nullptr; // 重建节点是当前节点，且subtree节点数量少于多线程重建阈值
  bool need_rebuild = allow_rebuild & Criterion_Check((*root));
  if (need_rebuild)
    Rebuild(root); // 主线程重建
  if ((*root) != nullptr)
    (*root)->working_flag = false;
  return;
}

template <typename PointType>
void KD_TREE<PointType>::Add_by_range(KD_TREE_NODE **root,
                                      BoxPointType boxpoint,
                                      bool allow_rebuild) {
  /**
   * \brief 类似delete_by_range
   */
  if ((*root) == nullptr)
    return;
  (*root)->working_flag = true;
  Push_Down(*root);
  if (boxpoint.vertex_max[0] <= (*root)->node_range_x[0] ||
      boxpoint.vertex_min[0] > (*root)->node_range_x[1])
    return;
  if (boxpoint.vertex_max[1] <= (*root)->node_range_y[0] ||
      boxpoint.vertex_min[1] > (*root)->node_range_y[1])
    return;
  if (boxpoint.vertex_max[2] <= (*root)->node_range_z[0] ||
      boxpoint.vertex_min[2] > (*root)->node_range_z[1])
    return;
  if (boxpoint.vertex_min[0] <= (*root)->node_range_x[0] &&
      boxpoint.vertex_max[0] > (*root)->node_range_x[1] &&
      boxpoint.vertex_min[1] <= (*root)->node_range_y[0] &&
      boxpoint.vertex_max[1] > (*root)->node_range_y[1] &&
      boxpoint.vertex_min[2] <= (*root)->node_range_z[0] &&
      boxpoint.vertex_max[2] > (*root)->node_range_z[1]) {
    (*root)->tree_deleted = false || (*root)->tree_downsample_deleted;
    (*root)->point_deleted = false || (*root)->point_downsample_deleted;
    (*root)->need_push_down_to_left = true;
    (*root)->need_push_down_to_right = true;
    (*root)->invalid_point_num = (*root)->down_del_num;
    return;
  }
  if (boxpoint.vertex_min[0] <= (*root)->point.x &&
      boxpoint.vertex_max[0] > (*root)->point.x &&
      boxpoint.vertex_min[1] <= (*root)->point.y &&
      boxpoint.vertex_max[1] > (*root)->point.y &&
      boxpoint.vertex_min[2] <= (*root)->point.z &&
      boxpoint.vertex_max[2] > (*root)->point.z) {
    (*root)->point_deleted = (*root)->point_downsample_deleted;
  }
  Operation_Logger_Type add_box_log;
  struct timespec Timeout;
  add_box_log.op = ADD_BOX;
  add_box_log.boxpoint = boxpoint;
  if ((Rebuild_Ptr == nullptr) || (*root)->left_son_ptr != *Rebuild_Ptr) {
    Add_by_range(&((*root)->left_son_ptr), boxpoint, allow_rebuild);
  } else {
    pthread_mutex_lock(&working_flag_mutex);
    Add_by_range(&((*root)->left_son_ptr), boxpoint, false);
    if (rebuild_flag) {
      pthread_mutex_lock(&rebuild_logger_mutex_lock);
      Rebuild_Logger.push(add_box_log);
      pthread_mutex_unlock(&rebuild_logger_mutex_lock);
    }
    pthread_mutex_unlock(&working_flag_mutex);
  }
  if ((Rebuild_Ptr == nullptr) || (*root)->right_son_ptr != *Rebuild_Ptr) {
    Add_by_range(&((*root)->right_son_ptr), boxpoint, allow_rebuild);
  } else {
    pthread_mutex_lock(&working_flag_mutex);
    Add_by_range(&((*root)->right_son_ptr), boxpoint, false);
    if (rebuild_flag) {
      pthread_mutex_lock(&rebuild_logger_mutex_lock);
      Rebuild_Logger.push(add_box_log);
      pthread_mutex_unlock(&rebuild_logger_mutex_lock);
    }
    pthread_mutex_unlock(&working_flag_mutex);
  }
  Update(*root);
  if (Rebuild_Ptr != nullptr && *Rebuild_Ptr == *root &&
      (*root)->TreeSize < Multi_Thread_Rebuild_Point_Num)
    Rebuild_Ptr = nullptr;
  bool need_rebuild = allow_rebuild & Criterion_Check((*root));
  if (need_rebuild)
    Rebuild(root);
  if ((*root) != nullptr)
    (*root)->working_flag = false;
  return;
}

template <typename PointType>
void KD_TREE<PointType>::Add_by_point(KD_TREE_NODE **root, PointType point,
                                      bool allow_rebuild, int father_axis) {
  if (*root == nullptr) {
    *root = new KD_TREE_NODE; // *root为指向根节点的指针
    InitTreeNode(*root);
    (*root)->point = point;
    (*root)->division_axis =
        (father_axis + 1) % 3; // 划分的轴为父节点的轴向后移一位
    Update(*root);
    return;
  } // Add_by_point最小单元
  (*root)->working_flag = true;
  Operation_Logger_Type add_log;
  struct timespec Timeout;
  add_log.op = ADD_POINT;
  add_log.point = point;
  Push_Down(*root);
  if (((*root)->division_axis == 0 && point.x < (*root)->point.x) ||
      ((*root)->division_axis == 1 && point.y < (*root)->point.y) ||
      ((*root)->division_axis == 2 &&
       point.z < (*root)->point.z)) { // 当该点处于左子树下时
    if ((Rebuild_Ptr == nullptr) || (*root)->left_son_ptr != *Rebuild_Ptr) {
      Add_by_point(&(*root)->left_son_ptr, point, allow_rebuild,
                   (*root)->division_axis);
    } else {
      pthread_mutex_lock(&working_flag_mutex);
      Add_by_point(&(*root)->left_son_ptr, point, false,
                   (*root)->division_axis); // allow_rebuild为false
      if (rebuild_flag) {
        pthread_mutex_lock(&rebuild_logger_mutex_lock);
        Rebuild_Logger.push(add_log);
        pthread_mutex_unlock(&rebuild_logger_mutex_lock);
      }
      pthread_mutex_unlock(&working_flag_mutex);
    }
  } else { // point落入右子树范围
    if ((Rebuild_Ptr == nullptr) || (*root)->right_son_ptr != *Rebuild_Ptr) {
      Add_by_point(&(*root)->right_son_ptr, point, allow_rebuild,
                   (*root)->division_axis);
    } else {
      pthread_mutex_lock(&working_flag_mutex);
      Add_by_point(&(*root)->right_son_ptr, point, false,
                   (*root)->division_axis);
      if (rebuild_flag) {
        pthread_mutex_lock(&rebuild_logger_mutex_lock);
        Rebuild_Logger.push(add_log);
        pthread_mutex_unlock(&rebuild_logger_mutex_lock);
      }
      pthread_mutex_unlock(&working_flag_mutex);
    }
  }
  Update(*root); // 更新当前节点的属性
  if (Rebuild_Ptr != nullptr && *Rebuild_Ptr == *root &&
      (*root)->TreeSize < Multi_Thread_Rebuild_Point_Num)
    Rebuild_Ptr = nullptr; // 当节点数量少于多线程重建阈值
  bool need_rebuild = allow_rebuild & Criterion_Check((*root));
  if (need_rebuild)
    Rebuild(root);
  if ((*root) != nullptr)
    (*root)->working_flag = false;
  return;
}

template <typename PointType>
void KD_TREE<PointType>::Search(KD_TREE_NODE *root, int k_nearest,
                                PointType point, MANUAL_HEAP &q,
                                double max_dist) {
  if (root == nullptr || root->tree_deleted)
    return;
  double cur_dist = calc_box_dist(root, point);
  double max_dist_sqr = max_dist * max_dist;
  if (cur_dist > max_dist_sqr) // max_dist_sqr为最大搜索距离
    return;
  int retval;
  if (root->need_push_down_to_left || root->need_push_down_to_right) {
    retval = pthread_mutex_trylock(&(root->push_down_mutex_lock)); // 返回0，表示加锁成功
    if (retval == 0) {
      Push_Down(root);
      pthread_mutex_unlock(&(root->push_down_mutex_lock));
    } else {
      pthread_mutex_lock(&(root->push_down_mutex_lock));
      pthread_mutex_unlock(&(root->push_down_mutex_lock));
    }
  }
  if (!root->point_deleted) { // 若当前节点没有被删除，则进行后续的比较
    float dist = calc_dist(point, root->point);
    if (dist <= max_dist_sqr && (q.size() < k_nearest || dist < q.top().dist)) { // 当前节点距离查询点距离小于优先级队列中堆顶的节点到查寻点距离
      if (q.size() >= k_nearest)
        q.pop(); // 首先将堆顶节点删除
      PointType_CMP current_point{root->point, dist};
      q.push(current_point);
    }
  }
  int cur_search_counter;
  float dist_left_node = calc_box_dist(root->left_son_ptr, point); // 进行递归运算point仍然为查询点
  float dist_right_node = calc_box_dist(root->right_son_ptr, point);
  if (q.size() < k_nearest ||
      dist_left_node < q.top().dist && dist_right_node < q.top().dist) { // 若查询到最近邻点个数小于设定阈值；或者查询点到节点左孩子和右孩子的距离都小于堆顶节点
    if (dist_left_node <= dist_right_node) { 
      // 左子树距离查询点更近，因此首先搜索左子树
      if (Rebuild_Ptr == nullptr || *Rebuild_Ptr != root->left_son_ptr) {
        Search(root->left_son_ptr, k_nearest, point, q, max_dist);
      } else {
        pthread_mutex_lock(&search_flag_mutex);
        while (search_mutex_counter == -1) {
          pthread_mutex_unlock(&search_flag_mutex);
          usleep(1);
          pthread_mutex_lock(&search_flag_mutex);
        }
        search_mutex_counter += 1;
        pthread_mutex_unlock(&search_flag_mutex);
        Search(root->left_son_ptr, k_nearest, point, q, max_dist);
        pthread_mutex_lock(&search_flag_mutex);
        search_mutex_counter -= 1;
        pthread_mutex_unlock(&search_flag_mutex);
      }
      // 左子树搜索完毕后，若搜索的近邻点数量仍不足且右子树比堆顶节点更近邻，则继续搜索右子树
      if (q.size() < k_nearest || dist_right_node < q.top().dist) {
        if (Rebuild_Ptr == nullptr || *Rebuild_Ptr != root->right_son_ptr) {
          Search(root->right_son_ptr, k_nearest, point, q, max_dist);
        } else {
          pthread_mutex_lock(&search_flag_mutex);
          while (search_mutex_counter == -1) {
            pthread_mutex_unlock(&search_flag_mutex);
            usleep(1);
            pthread_mutex_lock(&search_flag_mutex);
          }
          search_mutex_counter += 1;
          pthread_mutex_unlock(&search_flag_mutex);
          Search(root->right_son_ptr, k_nearest, point, q, max_dist);
          pthread_mutex_lock(&search_flag_mutex);
          search_mutex_counter -= 1;
          pthread_mutex_unlock(&search_flag_mutex);
        }
      }
    } else { // 右子树更靠近查询点
      // 先在右子树搜索
      if (Rebuild_Ptr == nullptr || *Rebuild_Ptr != root->right_son_ptr) {
        Search(root->right_son_ptr, k_nearest, point, q, max_dist);
      } else {
        pthread_mutex_lock(&search_flag_mutex);
        while (search_mutex_counter == -1) {
          pthread_mutex_unlock(&search_flag_mutex);
          usleep(1);
          pthread_mutex_lock(&search_flag_mutex);
        }
        search_mutex_counter += 1;
        pthread_mutex_unlock(&search_flag_mutex);
        Search(root->right_son_ptr, k_nearest, point, q, max_dist);
        pthread_mutex_lock(&search_flag_mutex);
        search_mutex_counter -= 1;
        pthread_mutex_unlock(&search_flag_mutex);
      }
      // 查询完右子树后，搜索到的最近邻点仍然不满足阈值，若此时左子树比堆顶节点更近邻，则可以查询左子树
      if (q.size() < k_nearest || dist_left_node < q.top().dist) {
        if (Rebuild_Ptr == nullptr || *Rebuild_Ptr != root->left_son_ptr) {
          Search(root->left_son_ptr, k_nearest, point, q, max_dist);
        } else {
          pthread_mutex_lock(&search_flag_mutex);
          while (search_mutex_counter == -1) {
            pthread_mutex_unlock(&search_flag_mutex);
            usleep(1);
            pthread_mutex_lock(&search_flag_mutex);
          }
          search_mutex_counter += 1;
          pthread_mutex_unlock(&search_flag_mutex);
          Search(root->left_son_ptr, k_nearest, point, q, max_dist);
          pthread_mutex_lock(&search_flag_mutex);
          search_mutex_counter -= 1;
          pthread_mutex_unlock(&search_flag_mutex);
        }
      }
    }
  } else { // 查询到的最近邻点数量大于阈值 || 左右孩子都比堆顶更靠近查询点
    if (dist_left_node < q.top().dist) { // 左孩子比堆顶更近邻
      if (Rebuild_Ptr == nullptr || *Rebuild_Ptr != root->left_son_ptr) {
        Search(root->left_son_ptr, k_nearest, point, q, max_dist);
      } else {
        pthread_mutex_lock(&search_flag_mutex);
        while (search_mutex_counter == -1) {
          pthread_mutex_unlock(&search_flag_mutex);
          usleep(1);
          pthread_mutex_lock(&search_flag_mutex);
        }
        search_mutex_counter += 1;
        pthread_mutex_unlock(&search_flag_mutex);
        Search(root->left_son_ptr, k_nearest, point, q, max_dist);
        pthread_mutex_lock(&search_flag_mutex);
        search_mutex_counter -= 1;
        pthread_mutex_unlock(&search_flag_mutex);
      }
    }
    if (dist_right_node < q.top().dist) { // 右孩子比堆顶更近邻
      if (Rebuild_Ptr == nullptr || *Rebuild_Ptr != root->right_son_ptr) {
        Search(root->right_son_ptr, k_nearest, point, q, max_dist);
      } else {
        pthread_mutex_lock(&search_flag_mutex);
        while (search_mutex_counter == -1) {
          pthread_mutex_unlock(&search_flag_mutex);
          usleep(1);
          pthread_mutex_lock(&search_flag_mutex);
        }
        search_mutex_counter += 1;
        pthread_mutex_unlock(&search_flag_mutex);
        Search(root->right_son_ptr, k_nearest, point, q, max_dist);
        pthread_mutex_lock(&search_flag_mutex);
        search_mutex_counter -= 1;
        pthread_mutex_unlock(&search_flag_mutex);
      }
    }
  }
  return;
}

template <typename PointType>
void KD_TREE<PointType>::Search_by_range(KD_TREE_NODE *root,
                                         BoxPointType boxpoint,
                                         PointVector &Storage) {
  /**
   * 1. 判断box和subtree在空间上是否存在交集
   * 2. 若subtree含于box，则直接采用flatten方法
   * 3. 若当前节点含于box，则push_back
   * 4. 递归左右子树，考虑多线程需要加锁
   */
  if (root == nullptr)
    return;
  Push_Down(root);
  if (boxpoint.vertex_max[0] <= root->node_range_x[0] ||
      boxpoint.vertex_min[0] > root->node_range_x[1])
    return;
  if (boxpoint.vertex_max[1] <= root->node_range_y[0] ||
      boxpoint.vertex_min[1] > root->node_range_y[1])
    return;
  if (boxpoint.vertex_max[2] <= root->node_range_z[0] ||
      boxpoint.vertex_min[2] > root->node_range_z[1])
    return;
  if (boxpoint.vertex_min[0] <= root->node_range_x[0] &&
      boxpoint.vertex_max[0] > root->node_range_x[1] &&
      boxpoint.vertex_min[1] <= root->node_range_y[0] &&
      boxpoint.vertex_max[1] > root->node_range_y[1] &&
      boxpoint.vertex_min[2] <= root->node_range_z[0] &&
      boxpoint.vertex_max[2] > root->node_range_z[1]) {
    flatten(root, Storage, NOT_RECORD);
    return;
  } // 若子树完全容于box中，则直接将root下所有节点保存
  if (boxpoint.vertex_min[0] <= root->point.x &&
      boxpoint.vertex_max[0] > root->point.x &&
      boxpoint.vertex_min[1] <= root->point.y &&
      boxpoint.vertex_max[1] > root->point.y &&
      boxpoint.vertex_min[2] <= root->point.z &&
      boxpoint.vertex_max[2] > root->point.z) {
    if (!root->point_deleted)
      Storage.push_back(root->point);
  } // 若当前节点坐标含于box且有效
  if ((Rebuild_Ptr == nullptr) || root->left_son_ptr != *Rebuild_Ptr) {
    Search_by_range(root->left_son_ptr, boxpoint, Storage);
  } else {
    pthread_mutex_lock(&search_flag_mutex);
    Search_by_range(root->left_son_ptr, boxpoint, Storage);
    pthread_mutex_unlock(&search_flag_mutex);
  }
  if ((Rebuild_Ptr == nullptr) || root->right_son_ptr != *Rebuild_Ptr) {
    Search_by_range(root->right_son_ptr, boxpoint, Storage);
  } else {
    pthread_mutex_lock(&search_flag_mutex);
    Search_by_range(root->right_son_ptr, boxpoint, Storage);
    pthread_mutex_unlock(&search_flag_mutex);
  }
  return;
}

template <typename PointType>
void KD_TREE<PointType>::Search_by_radius(KD_TREE_NODE *root, PointType point,
                                          float radius, PointVector &Storage) {
  if (root == nullptr)
    return;
  Push_Down(root);
  PointType range_center; // range_box的中心点
  range_center.x = (root->node_range_x[0] + root->node_range_x[1]) * 0.5;
  range_center.y = (root->node_range_y[0] + root->node_range_y[1]) * 0.5;
  range_center.z = (root->node_range_z[0] + root->node_range_z[1]) * 0.5;
  float dist = sqrt(calc_dist(range_center, point)); // 计算point与当前节点的range_box中心点的距离
  if (dist > radius + sqrt(root->radius_sq)) // 查询点的查询半径和当前节点下的subtree在空间上没有交集
    return;
  if (dist <= radius - sqrt(root->radius_sq)) { // 当前节点的subtree的range范围含于查询点查询半径范围（内切于query ball）
    flatten(root, Storage, NOT_RECORD);
    return;
  }
  if (!root->point_deleted &&
      calc_dist(root->point, point) <= radius * radius) { // 计算当前节点距离查询点距离
    Storage.push_back(root->point); // 当前节点在查询半径内
  }
  // 搜索左子树
  if ((Rebuild_Ptr == nullptr) || root->left_son_ptr != *Rebuild_Ptr) {
    Search_by_radius(root->left_son_ptr, point, radius, Storage);
  } else {
    pthread_mutex_lock(&search_flag_mutex);
    Search_by_radius(root->left_son_ptr, point, radius, Storage);
    pthread_mutex_unlock(&search_flag_mutex);
  }
  // 搜索右子树
  if ((Rebuild_Ptr == nullptr) || root->right_son_ptr != *Rebuild_Ptr) {
    Search_by_radius(root->right_son_ptr, point, radius, Storage);
  } else {
    pthread_mutex_lock(&search_flag_mutex);
    Search_by_radius(root->right_son_ptr, point, radius, Storage);
    pthread_mutex_unlock(&search_flag_mutex);
  }
  return;
}

template <typename PointType>
bool KD_TREE<PointType>::Criterion_Check(
    KD_TREE_NODE *root) { // 检查平衡标准,bal和del
  if (root->TreeSize <= Minimal_Unbalanced_Tree_Size) {
    return false;
  }
  float balance_evaluation = 0.0f;
  float delete_evaluation = 0.0f;
  KD_TREE_NODE *son_ptr = root->left_son_ptr;
  if (son_ptr == nullptr)
    son_ptr = root->right_son_ptr;
  delete_evaluation = float(root->invalid_point_num) / root->TreeSize;
  balance_evaluation = float(son_ptr->TreeSize) / (root->TreeSize - 1);
  if (delete_evaluation > delete_criterion_param) {
    return true;
  } // 无效点超过阈值
  if (balance_evaluation > balance_criterion_param ||
      balance_evaluation < 1 - balance_criterion_param) {
    return true;
  } // 左右子树不平衡，数量差距过大
  return false;
}

template <typename PointType>
void KD_TREE<PointType>::Push_Down(KD_TREE_NODE *root) {
  /**
   * 1. 对于box-wise
   * updates，在进行delete操作时，需要递归的更新deleted和treedeleted属性，常规操作不够高效
   * 2. 于是采用Push_Down这种further lazy
   * strategy(对应论文中)，kdtree每个节点也增加了pushdown这个属性
   * 3.
   * Push_Down操作copy当前节点T的deleted、treedeleted以及pushdown属性到其sons（仅限于下一层的sons,
   * not further offsprings）
   */
  if (root == nullptr)
    return;
  Operation_Logger_Type operation;
  operation.op = PUSH_DOWN;
  operation.tree_deleted = root->tree_deleted; // 当前节点的tree_deleted
  operation.tree_downsample_deleted =
      root->tree_downsample_deleted; // tree_downsample_deleted表示downsample删除的节点
  if (root->need_push_down_to_left &&
      root->left_son_ptr != nullptr) { // 若左子树需要下推且存在左子树
    if (Rebuild_Ptr == nullptr ||
        *Rebuild_Ptr !=
            root->left_son_ptr) { // 若重建的节点为空或者不为当前的左子节点
      root->left_son_ptr->tree_downsample_deleted |=
          root->tree_downsample_deleted; // 位或操作
      root->left_son_ptr->point_downsample_deleted |=
          root->tree_downsample_deleted;
      root->left_son_ptr->tree_deleted =
          root->tree_deleted || root->left_son_ptr->tree_downsample_deleted;
      root->left_son_ptr->point_deleted =
          root->left_son_ptr->tree_deleted ||
          root->left_son_ptr->point_downsample_deleted;
      if (root->tree_downsample_deleted)
        root->left_son_ptr->down_del_num =
            root->left_son_ptr
                ->TreeSize; // 若root节点下采样时被删除，则root节点的子节点down_del_num置为左子树的TreeSize
      if (root->tree_deleted)
        root->left_son_ptr->invalid_point_num = root->left_son_ptr->TreeSize;
      else
        root->left_son_ptr->invalid_point_num =
            root->left_son_ptr
                ->down_del_num; // invalid_point_num等于下采样被删除的点个数
      root->left_son_ptr->need_push_down_to_left = true;
      root->left_son_ptr->need_push_down_to_right = true;
      root->need_push_down_to_left =
          false; // push_down结束后当前节点的pushdown属性置false
    } else {     // 若重建的是当前节点
      pthread_mutex_lock(
          &working_flag_mutex); // 若当前节点在重建中，则需要线程加锁
      root->left_son_ptr->tree_downsample_deleted |=
          root->tree_downsample_deleted;
      root->left_son_ptr->point_downsample_deleted |=
          root->tree_downsample_deleted;
      root->left_son_ptr->tree_deleted =
          root->tree_deleted || root->left_son_ptr->tree_downsample_deleted;
      root->left_son_ptr->point_deleted =
          root->left_son_ptr->tree_deleted ||
          root->left_son_ptr->point_downsample_deleted;
      if (root->tree_downsample_deleted)
        root->left_son_ptr->down_del_num = root->left_son_ptr->TreeSize;
      if (root->tree_deleted)
        root->left_son_ptr->invalid_point_num = root->left_son_ptr->TreeSize;
      else
        root->left_son_ptr->invalid_point_num =
            root->left_son_ptr->down_del_num;
      root->left_son_ptr->need_push_down_to_left = true;
      root->left_son_ptr->need_push_down_to_right = true;
      if (rebuild_flag) { // 若rebuild_flag为true，暂存下面的操作
        pthread_mutex_lock(&rebuild_logger_mutex_lock);
        Rebuild_Logger.push(operation);
        pthread_mutex_unlock(&rebuild_logger_mutex_lock);
      }
      root->need_push_down_to_left = false;
      pthread_mutex_unlock(&working_flag_mutex);
    }
  }
  if (root->need_push_down_to_right &&
      root->right_son_ptr != nullptr) { // 右子树需要push_down且存在右子树
    if (Rebuild_Ptr == nullptr || *Rebuild_Ptr != root->right_son_ptr) {
      root->right_son_ptr->tree_downsample_deleted |=
          root->tree_downsample_deleted;
      root->right_son_ptr->point_downsample_deleted |=
          root->tree_downsample_deleted;
      root->right_son_ptr->tree_deleted =
          root->tree_deleted || root->right_son_ptr->tree_downsample_deleted;
      root->right_son_ptr->point_deleted =
          root->right_son_ptr->tree_deleted ||
          root->right_son_ptr->point_downsample_deleted;
      if (root->tree_downsample_deleted)
        root->right_son_ptr->down_del_num = root->right_son_ptr->TreeSize;
      if (root->tree_deleted)
        root->right_son_ptr->invalid_point_num = root->right_son_ptr->TreeSize;
      else
        root->right_son_ptr->invalid_point_num =
            root->right_son_ptr->down_del_num;
      root->right_son_ptr->need_push_down_to_left = true;
      root->right_son_ptr->need_push_down_to_right = true;
      root->need_push_down_to_right = false;
    } else {
      pthread_mutex_lock(&working_flag_mutex);
      root->right_son_ptr->tree_downsample_deleted |=
          root->tree_downsample_deleted;
      root->right_son_ptr->point_downsample_deleted |=
          root->tree_downsample_deleted;
      root->right_son_ptr->tree_deleted =
          root->tree_deleted || root->right_son_ptr->tree_downsample_deleted;
      root->right_son_ptr->point_deleted =
          root->right_son_ptr->tree_deleted ||
          root->right_son_ptr->point_downsample_deleted;
      if (root->tree_downsample_deleted)
        root->right_son_ptr->down_del_num = root->right_son_ptr->TreeSize;
      if (root->tree_deleted)
        root->right_son_ptr->invalid_point_num = root->right_son_ptr->TreeSize;
      else
        root->right_son_ptr->invalid_point_num =
            root->right_son_ptr->down_del_num;
      root->right_son_ptr->need_push_down_to_left = true;
      root->right_son_ptr->need_push_down_to_right = true;
      if (rebuild_flag) { // 若正在重建，则暂存操作
        pthread_mutex_lock(&rebuild_logger_mutex_lock);
        Rebuild_Logger.push(operation);
        pthread_mutex_unlock(&rebuild_logger_mutex_lock);
      }
      root->need_push_down_to_right = false;
      pthread_mutex_unlock(&working_flag_mutex);
    }
  }
  return;
}

template <typename PointType>
void KD_TREE<PointType>::Update(
    KD_TREE_NODE *
        root) { // 更新重建的subtree相关属性（treesize, invalidnum, range etc.）
  /**
   * 1. 对应论文中的Pullup操作
   * 2.
   * 总结位于当前节点T下subtree的相关属性，包括treesize、invalidnum、down_del_num、range
   * 3. 也会更新tree_downsample_deleted、tree_deleted。
   * 4.
   * 感觉是与push_down对应，psuh_down向下更新，而pullup是向上，将subtree的信息总结到但前节点
   */
  KD_TREE_NODE *left_son_ptr = root->left_son_ptr; // 左孩子
  KD_TREE_NODE *right_son_ptr = root->right_son_ptr;
  float tmp_range_x[2] = {INFINITY, -INFINITY};
  float tmp_range_y[2] = {INFINITY, -INFINITY};
  float tmp_range_z[2] = {INFINITY, -INFINITY};
  // Update Tree Size
  if (left_son_ptr != nullptr && right_son_ptr != nullptr) { // 左、右孩子都存在
    root->TreeSize = left_son_ptr->TreeSize + right_son_ptr->TreeSize + 1;
    root->invalid_point_num = left_son_ptr->invalid_point_num +
                              right_son_ptr->invalid_point_num +
                              (root->point_deleted ? 1 : 0);
    root->down_del_num = left_son_ptr->down_del_num +
                         right_son_ptr->down_del_num +
                         (root->point_downsample_deleted ? 1 : 0);
    root->tree_downsample_deleted = left_son_ptr->tree_downsample_deleted &
                                    right_son_ptr->tree_downsample_deleted &
                                    root->point_downsample_deleted;
    root->tree_deleted = left_son_ptr->tree_deleted &&
                         right_son_ptr->tree_deleted && root->point_deleted;
    if (root->tree_deleted ||
        (!left_son_ptr->tree_deleted && !right_son_ptr->tree_deleted &&
         !root->point_deleted)) { // 当该节点的subtree被删除或者左右子树和当前节点都未被删除，只需要执行下面的逻辑
      tmp_range_x[0] = min(
          min(left_son_ptr->node_range_x[0], right_son_ptr->node_range_x[0]),
          root->point.x);
      tmp_range_x[1] = max(
          max(left_son_ptr->node_range_x[1], right_son_ptr->node_range_x[1]),
          root->point.x);
      tmp_range_y[0] = min(
          min(left_son_ptr->node_range_y[0], right_son_ptr->node_range_y[0]),
          root->point.y);
      tmp_range_y[1] = max(
          max(left_son_ptr->node_range_y[1], right_son_ptr->node_range_y[1]),
          root->point.y);
      tmp_range_z[0] = min(
          min(left_son_ptr->node_range_z[0], right_son_ptr->node_range_z[0]),
          root->point.z);
      tmp_range_z[1] = max(
          max(left_son_ptr->node_range_z[1], right_son_ptr->node_range_z[1]),
          root->point.z);
    } else { // root->tree_deleted=false
             // left_son_ptr->tree_deleted、right_son_ptr->tree_deleted、root->point_deleted中有被删除的
      if (!left_son_ptr->tree_deleted) { // 左子树未被删除
        tmp_range_x[0] = min(tmp_range_x[0], left_son_ptr->node_range_x[0]);
        tmp_range_x[1] = max(tmp_range_x[1], left_son_ptr->node_range_x[1]);
        tmp_range_y[0] = min(tmp_range_y[0], left_son_ptr->node_range_y[0]);
        tmp_range_y[1] = max(tmp_range_y[1], left_son_ptr->node_range_y[1]);
        tmp_range_z[0] = min(tmp_range_z[0], left_son_ptr->node_range_z[0]);
        tmp_range_z[1] = max(tmp_range_z[1], left_son_ptr->node_range_z[1]);
      }
      if (!right_son_ptr->tree_deleted) { // 右子树未被删除
        tmp_range_x[0] = min(tmp_range_x[0], right_son_ptr->node_range_x[0]);
        tmp_range_x[1] = max(tmp_range_x[1], right_son_ptr->node_range_x[1]);
        tmp_range_y[0] = min(tmp_range_y[0], right_son_ptr->node_range_y[0]);
        tmp_range_y[1] = max(tmp_range_y[1], right_son_ptr->node_range_y[1]);
        tmp_range_z[0] = min(tmp_range_z[0], right_son_ptr->node_range_z[0]);
        tmp_range_z[1] = max(tmp_range_z[1], right_son_ptr->node_range_z[1]);
      }
      if (!root->point_deleted) { // 当前节点未被删除
        tmp_range_x[0] = min(tmp_range_x[0], root->point.x);
        tmp_range_x[1] = max(tmp_range_x[1], root->point.x);
        tmp_range_y[0] = min(tmp_range_y[0], root->point.y);
        tmp_range_y[1] = max(tmp_range_y[1], root->point.y);
        tmp_range_z[0] = min(tmp_range_z[0], root->point.z);
        tmp_range_z[1] = max(tmp_range_z[1], root->point.z);
      }
      // 上述操作可以依次操作，可以更新tmp_range
    }
  } else if (left_son_ptr != nullptr) { // 只有左孩子存在
    root->TreeSize = left_son_ptr->TreeSize + 1;
    root->invalid_point_num =
        left_son_ptr->invalid_point_num + (root->point_deleted ? 1 : 0);
    root->down_del_num =
        left_son_ptr->down_del_num + (root->point_downsample_deleted ? 1 : 0);
    root->tree_downsample_deleted =
        left_son_ptr->tree_downsample_deleted & root->point_downsample_deleted;
    root->tree_deleted = left_son_ptr->tree_deleted && root->point_deleted;
    if (root->tree_deleted ||
        (!left_son_ptr->tree_deleted && !root->point_deleted)) {
      tmp_range_x[0] = min(left_son_ptr->node_range_x[0], root->point.x);
      tmp_range_x[1] = max(left_son_ptr->node_range_x[1], root->point.x);
      tmp_range_y[0] = min(left_son_ptr->node_range_y[0], root->point.y);
      tmp_range_y[1] = max(left_son_ptr->node_range_y[1], root->point.y);
      tmp_range_z[0] = min(left_son_ptr->node_range_z[0], root->point.z);
      tmp_range_z[1] = max(left_son_ptr->node_range_z[1], root->point.z);
    } else {
      if (!left_son_ptr->tree_deleted) {
        tmp_range_x[0] = min(tmp_range_x[0], left_son_ptr->node_range_x[0]);
        tmp_range_x[1] = max(tmp_range_x[1], left_son_ptr->node_range_x[1]);
        tmp_range_y[0] = min(tmp_range_y[0], left_son_ptr->node_range_y[0]);
        tmp_range_y[1] = max(tmp_range_y[1], left_son_ptr->node_range_y[1]);
        tmp_range_z[0] = min(tmp_range_z[0], left_son_ptr->node_range_z[0]);
        tmp_range_z[1] = max(tmp_range_z[1], left_son_ptr->node_range_z[1]);
      }
      if (!root->point_deleted) {
        tmp_range_x[0] = min(tmp_range_x[0], root->point.x);
        tmp_range_x[1] = max(tmp_range_x[1], root->point.x);
        tmp_range_y[0] = min(tmp_range_y[0], root->point.y);
        tmp_range_y[1] = max(tmp_range_y[1], root->point.y);
        tmp_range_z[0] = min(tmp_range_z[0], root->point.z);
        tmp_range_z[1] = max(tmp_range_z[1], root->point.z);
      }
    }

  } else if (right_son_ptr != nullptr) { // 只有右孩子存在
    root->TreeSize = right_son_ptr->TreeSize + 1;
    root->invalid_point_num =
        right_son_ptr->invalid_point_num + (root->point_deleted ? 1 : 0);
    root->down_del_num =
        right_son_ptr->down_del_num + (root->point_downsample_deleted ? 1 : 0);
    root->tree_downsample_deleted =
        right_son_ptr->tree_downsample_deleted & root->point_downsample_deleted;
    root->tree_deleted = right_son_ptr->tree_deleted && root->point_deleted;
    if (root->tree_deleted ||
        (!right_son_ptr->tree_deleted && !root->point_deleted)) {
      tmp_range_x[0] = min(right_son_ptr->node_range_x[0], root->point.x);
      tmp_range_x[1] = max(right_son_ptr->node_range_x[1], root->point.x);
      tmp_range_y[0] = min(right_son_ptr->node_range_y[0], root->point.y);
      tmp_range_y[1] = max(right_son_ptr->node_range_y[1], root->point.y);
      tmp_range_z[0] = min(right_son_ptr->node_range_z[0], root->point.z);
      tmp_range_z[1] = max(right_son_ptr->node_range_z[1], root->point.z);
    } else {
      if (!right_son_ptr->tree_deleted) {
        tmp_range_x[0] = min(tmp_range_x[0], right_son_ptr->node_range_x[0]);
        tmp_range_x[1] = max(tmp_range_x[1], right_son_ptr->node_range_x[1]);
        tmp_range_y[0] = min(tmp_range_y[0], right_son_ptr->node_range_y[0]);
        tmp_range_y[1] = max(tmp_range_y[1], right_son_ptr->node_range_y[1]);
        tmp_range_z[0] = min(tmp_range_z[0], right_son_ptr->node_range_z[0]);
        tmp_range_z[1] = max(tmp_range_z[1], right_son_ptr->node_range_z[1]);
      }
      if (!root->point_deleted) {
        tmp_range_x[0] = min(tmp_range_x[0], root->point.x);
        tmp_range_x[1] = max(tmp_range_x[1], root->point.x);
        tmp_range_y[0] = min(tmp_range_y[0], root->point.y);
        tmp_range_y[1] = max(tmp_range_y[1], root->point.y);
        tmp_range_z[0] = min(tmp_range_z[0], root->point.z);
        tmp_range_z[1] = max(tmp_range_z[1], root->point.z);
      }
    }
  } else { // 没有孩子
    root->TreeSize = 1;
    root->invalid_point_num = (root->point_deleted ? 1 : 0);
    root->down_del_num = (root->point_downsample_deleted ? 1 : 0);
    root->tree_downsample_deleted = root->point_downsample_deleted;
    root->tree_deleted = root->point_deleted;
    tmp_range_x[0] = root->point.x;
    tmp_range_x[1] = root->point.x;
    tmp_range_y[0] = root->point.y;
    tmp_range_y[1] = root->point.y;
    tmp_range_z[0] = root->point.z;
    tmp_range_z[1] = root->point.z;
  }
  memcpy(root->node_range_x, tmp_range_x, sizeof(tmp_range_x));
  memcpy(root->node_range_y, tmp_range_y, sizeof(tmp_range_y));
  memcpy(root->node_range_z, tmp_range_z, sizeof(tmp_range_z));
  float x_L = (root->node_range_x[1] - root->node_range_x[0]) * 0.5;
  float y_L = (root->node_range_y[1] - root->node_range_y[0]) * 0.5;
  float z_L = (root->node_range_z[1] - root->node_range_z[0]) * 0.5;
  root->radius_sq =
      x_L * x_L + y_L * y_L +
      z_L * z_L; // 半径的平方，box的最近最远点的连线距离的一半的平方
  if (left_son_ptr != nullptr)
    left_son_ptr->father_ptr = root; // 更新father_ptr
  if (right_son_ptr != nullptr)
    right_son_ptr->father_ptr = root;
  if (root == Root_Node &&
      root->TreeSize > 3) { // 当该节点为根节点，且treesize>3时;
                            // 统计平衡因子，alpha_bal和alpha_del
    KD_TREE_NODE *son_ptr = root->left_son_ptr;
    if (son_ptr == nullptr)
      son_ptr = root->right_son_ptr;
    float tmp_bal =
        float(son_ptr->TreeSize) /
        (root->TreeSize - 1); // 计算平衡因子，左/右子树节点数除以子树节点数之和
    root->alpha_del =
        float(root->invalid_point_num) /
        root->TreeSize; // 计算删除因子，删除点数除以节点数,取值范围为(0~1)
    root->alpha_bal = (tmp_bal >= 0.5 - EPSS)
                          ? tmp_bal
                          : 1 - tmp_bal; // alpha_bal取值范围为(0.5~1)
  }
  return;
}

template <typename PointType>
void KD_TREE<PointType>::flatten(KD_TREE_NODE *root, PointVector &Storage,
                                 delete_point_storage_set storage_type) {
  if (root == nullptr)
    return; // 终止递归
  Push_Down(root);
  if (!root->point_deleted) { // 若为有效节点，则将该节点添加到Storage
    Storage.push_back(root->point);
  }
  flatten(root->left_son_ptr, Storage, storage_type);  // 遍历左子树
  flatten(root->right_son_ptr, Storage, storage_type); // 遍历右子树
  switch (storage_type) {
  case NOT_RECORD:
    break;
  case DELETE_POINTS_REC:
    if (root->point_deleted && !root->point_downsample_deleted) {
      Points_deleted.push_back(root->point); // 添加入点删除列表
    }
    break;
  case MULTI_THREAD_REC:
    if (root->point_deleted && !root->point_downsample_deleted) {
      Multithread_Points_deleted.push_back(root->point);
    }
    break;
  default:
    break;
  }
  return;
}

template <typename PointType>
void KD_TREE<PointType>::delete_tree_nodes(KD_TREE_NODE **root) {
  if (*root == nullptr)
    return;
  Push_Down(*root); // 将当前节点的一些属性copy到子节点
  delete_tree_nodes(&(*root)->left_son_ptr);
  delete_tree_nodes(&(*root)->right_son_ptr);

  pthread_mutex_destroy(&(*root)->push_down_mutex_lock);
  delete *root;
  *root = nullptr;

  return;
}

template <typename PointType>
bool KD_TREE<PointType>::same_point(PointType a,
                                    PointType b) { // 判断两个点是否相同
  return (fabs(a.x - b.x) < EPSS && fabs(a.y - b.y) < EPSS &&
          fabs(a.z - b.z) < EPSS);
}

template <typename PointType>
float KD_TREE<PointType>::calc_dist(PointType a,
                                    PointType b) { // 计算点之间的平方距离
  float dist = 0.0f;
  dist = (a.x - b.x) * (a.x - b.x) + (a.y - b.y) * (a.y - b.y) +
         (a.z - b.z) * (a.z - b.z);
  return dist;
}

template <typename PointType>
float KD_TREE<PointType>::calc_box_dist(
    KD_TREE_NODE *node, PointType point) { 
  // 计算当前点与节点range_box的距离
  /**
   * 1. 若point落入node的range_box，则min_dist为0
   * 2. 计算point与range_box的最大最远点距离
   */
  if (node == nullptr)
    return INFINITY;
  float min_dist = 0.0;
  if (point.x < node->node_range_x[0])
    min_dist +=
        (point.x - node->node_range_x[0]) * (point.x - node->node_range_x[0]); // point.x小于range_box.min_x
  if (point.x > node->node_range_x[1])
    min_dist +=
        (point.x - node->node_range_x[1]) * (point.x - node->node_range_x[1]); // point.x大于range_box.max_x
  if (point.y < node->node_range_y[0])
    min_dist +=
        (point.y - node->node_range_y[0]) * (point.y - node->node_range_y[0]); // point.y小于range_box.min_y
  if (point.y > node->node_range_y[1])
    min_dist +=
        (point.y - node->node_range_y[1]) * (point.y - node->node_range_y[1]); // point.y小于range_box.max_x
  if (point.z < node->node_range_z[0])
    min_dist +=
        (point.z - node->node_range_z[0]) * (point.z - node->node_range_z[0]); // point.z小于range_box.min_z
  if (point.z > node->node_range_z[1])
    min_dist +=
        (point.z - node->node_range_z[1]) * (point.z - node->node_range_z[1]); // point.z小于range_box.max_x
  return min_dist;
}

template <typename PointType>
bool KD_TREE<PointType>::point_cmp_x(PointType a, PointType b) {
  return a.x < b.x;
}
template <typename PointType>
bool KD_TREE<PointType>::point_cmp_y(PointType a, PointType b) {
  return a.y < b.y;
}
template <typename PointType>
bool KD_TREE<PointType>::point_cmp_z(PointType a, PointType b) {
  return a.z < b.z;
}

// manual queue
template <typename T> void MANUAL_Q<T>::clear() {
  head = 0;
  tail = 0;
  counter = 0;
  is_empty = true;
  return;
} // 头、尾、计数均置0

template <typename T> void MANUAL_Q<T>::pop() {
  if (counter == 0)
    return;
  head++;
  head %= Q_LEN;
  counter--;
  if (counter == 0)
    is_empty = true;
  return;
} // 头指针后移

template <typename T> T MANUAL_Q<T>::front() { return q[head]; }

template <typename T> T MANUAL_Q<T>::back() { return q[tail]; }

template <typename T> void MANUAL_Q<T>::push(T op) {
  q[tail] = op;
  counter++;
  if (is_empty)
    is_empty = false;
  tail++;
  tail %= Q_LEN;
} // 尾指针后移，在尾部插入元素

template <typename T> bool MANUAL_Q<T>::empty() { return is_empty; }

template <typename T> int MANUAL_Q<T>::size() { return counter; }

template class KD_TREE<ikdTree_PointType>;
template class KD_TREE<pcl::PointXYZ>;
template class KD_TREE<pcl::PointXYZI>;
template class KD_TREE<pcl::PointXYZINormal>;