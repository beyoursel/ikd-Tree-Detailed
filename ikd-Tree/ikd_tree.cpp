#include "ikd_tree.h"

/*
Description: ikd-Tree: an incremental k-d tree for robotic applications
Author: Yixi Cai
email: yixicai@connect.hku.hk
*/

template <typename PointType>
KD_TREE<PointType>::KD_TREE(float delete_param, float balance_param,
                            float box_length) {
  delete_criterion_param_ = delete_param;
  balance_criterion_param_ = balance_param;
  downsample_size_ = box_length;
  root_node_ = nullptr;
  static_root_node_ = nullptr;
}

template <typename PointType>
KD_TREE<PointType>::~KD_TREE() {
  DeleteTreeNodes(&root_node_);
  PointVector().swap(pcl_storage_);
}

template <typename PointType>
void KD_TREE<PointType>::SetDeleteCriterionParam(float delete_param) {
  delete_criterion_param_ = delete_param;
}

template <typename PointType>
void KD_TREE<PointType>::SetBalanceCriterionParam(float balance_param) {
  balance_criterion_param_ = balance_param;
}

template <typename PointType>
void KD_TREE<PointType>::SetDownsampleParam(float downsample_param) {
  downsample_size_ = downsample_param;
}

template <typename PointType>
void KD_TREE<PointType>::InitializeKDTree(float delete_param,
                                          float balance_param,
                                          float box_length) {
  SetDeleteCriterionParam(delete_param);
  SetBalanceCriterionParam(balance_param);
  SetDownsampleParam(box_length);
}

template <typename PointType>
void KD_TREE<PointType>::InitTreeNode(KD_TREE_NODE* root) {
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
  root->tree_size = 0;
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

template <typename PointType>
int KD_TREE<PointType>::Size() {
  int s = 0;
  if (root_node_ != nullptr) {
    return root_node_->tree_size;
  } else {
    return 0;
  }
}

template <typename PointType>
BoxPointType KD_TREE<PointType>::TreeRange() {
  BoxPointType range;
  if (root_node_ != nullptr) {
    range.vertex_min[0] = root_node_->node_range_x[0];
    range.vertex_min[1] = root_node_->node_range_y[0];
    range.vertex_min[2] = root_node_->node_range_z[0];
    range.vertex_max[0] = root_node_->node_range_x[1];
    range.vertex_max[1] = root_node_->node_range_y[1];
    range.vertex_max[2] = root_node_->node_range_z[1];
  } else {
    memset(&range, 0, sizeof(range));
  }
  return range;
}

template <typename PointType>
int KD_TREE<PointType>::ValidNum() {
  int s = 0;
  if (root_node_ != nullptr)
    return (root_node_->tree_size - root_node_->invalid_point_num);
  else
    return 0;
}

template <typename PointType>
void KD_TREE<PointType>::RootAlpha(float& alpha_bal, float& alpha_del) {
  alpha_bal = root_node_->alpha_bal;
  alpha_del = root_node_->alpha_del;
}

template <typename PointType>
void KD_TREE<PointType>::Build(PointVector point_cloud) {
  if (root_node_ != nullptr) {
    DeleteTreeNodes(&root_node_);
  }
  if (point_cloud.size() == 0) return;
  static_root_node_ = new KD_TREE_NODE;
  InitTreeNode(static_root_node_);
  BuildTree(&static_root_node_->left_son_ptr, 0, point_cloud.size() - 1,
            point_cloud);
  Update(static_root_node_);
  static_root_node_->tree_size = 0;
  root_node_ = static_root_node_->left_son_ptr;
}

template <typename PointType>
void KD_TREE<PointType>::NearestSearch(PointType point, int k_nearest,
                                        PointVector& Nearest_Points,
                                        std::vector<float>& Point_Distance,
                                        double max_dist) {
  MANUAL_HEAP q(2 * k_nearest);
  q.clear();
  std::vector<float>().swap(Point_Distance);
  Search(root_node_, k_nearest, point, q, max_dist);
  int k_found = std::min(k_nearest, int(q.size()));
  PointVector().swap(Nearest_Points);
  std::vector<float>().swap(Point_Distance);
  for (int i = 0; i < k_found; i++) {
    Nearest_Points.insert(Nearest_Points.begin(), q.top().point);
    Point_Distance.insert(Point_Distance.begin(), q.top().dist);
    q.pop();
  }
  return;
}

template <typename PointType>
void KD_TREE<PointType>::BoxSearch(const BoxPointType& Box_of_Point,
                                    PointVector& Storage) {
  Storage.clear();
  SearchByRange(root_node_, Box_of_Point, Storage);
}

template <typename PointType>
void KD_TREE<PointType>::RadiusSearch(PointType point, const float radius,
                                       PointVector& Storage) {
  Storage.clear();
  SearchByRadius(root_node_, point, radius, Storage);
}

template <typename PointType>
int KD_TREE<PointType>::AddPoints(PointVector& PointToAdd,
                                   bool downsample_on) {
  int NewPointSize = PointToAdd.size();
  int tree_size = Size();
  BoxPointType Box_of_Point;
  PointType downsample_result, mid_point;
  bool downsample_switch = downsample_on && DOWNSAMPLE_SWITCH;
  float min_dist, tmp_dist;
  int tmp_counter = 0;
  for (int i = 0; i < PointToAdd.size(); i++) {
    if (downsample_switch) {
      Box_of_Point.vertex_min[0] =
          floor(PointToAdd[i].x / downsample_size_) * downsample_size_;
      Box_of_Point.vertex_max[0] = Box_of_Point.vertex_min[0] + downsample_size_;
      Box_of_Point.vertex_min[1] =
          floor(PointToAdd[i].y / downsample_size_) * downsample_size_;
      Box_of_Point.vertex_max[1] = Box_of_Point.vertex_min[1] + downsample_size_;
      Box_of_Point.vertex_min[2] =
          floor(PointToAdd[i].z / downsample_size_) * downsample_size_;
      Box_of_Point.vertex_max[2] = Box_of_Point.vertex_min[2] + downsample_size_;
      mid_point.x =
          Box_of_Point.vertex_min[0] +
          (Box_of_Point.vertex_max[0] - Box_of_Point.vertex_min[0]) / 2.0;
      mid_point.y =
          Box_of_Point.vertex_min[1] +
          (Box_of_Point.vertex_max[1] - Box_of_Point.vertex_min[1]) / 2.0;
      mid_point.z =
          Box_of_Point.vertex_min[2] +
          (Box_of_Point.vertex_max[2] - Box_of_Point.vertex_min[2]) / 2.0;
      PointVector().swap(downsample_storage_);
      SearchByRange(root_node_, Box_of_Point, downsample_storage_);
      min_dist = CalcDist(PointToAdd[i], mid_point);
      downsample_result = PointToAdd[i];
      for (int index = 0; index < downsample_storage_.size(); index++) {
        tmp_dist = CalcDist(downsample_storage_[index], mid_point);
        if (tmp_dist < min_dist) {
          min_dist = tmp_dist;
          downsample_result = downsample_storage_[index];
        }
      }
      if (downsample_storage_.size() > 1 ||
          SamePoint(PointToAdd[i], downsample_result)) {
        if (downsample_storage_.size() > 0)
          DeleteByRange(&root_node_, Box_of_Point, true, true);
        AddByPoint(&root_node_, downsample_result, true,
                     root_node_->division_axis);
        tmp_counter++;
      }
    } else {
      AddByPoint(&root_node_, PointToAdd[i], true, root_node_->division_axis);
    }
  }
  return tmp_counter;
}

template <typename PointType>
int KD_TREE<PointType>::DeleteByRange(KD_TREE_NODE **root, BoxPointType boxpoint, bool allow_rebuild, bool is_downsample)
{
    if ((*root) == nullptr || (*root)->tree_deleted)
        return 0;
    (*root)->working_flag = true;
    PushDown(*root);
    int tmp_counter = 0;
    if (boxpoint.vertex_max[0] <= (*root)->node_range_x[0] || boxpoint.vertex_min[0] > (*root)->node_range_x[1])
        return 0;
    if (boxpoint.vertex_max[1] <= (*root)->node_range_y[0] || boxpoint.vertex_min[1] > (*root)->node_range_y[1])
        return 0;
    if (boxpoint.vertex_max[2] <= (*root)->node_range_z[0] || boxpoint.vertex_min[2] > (*root)->node_range_z[1])
        return 0;
    if (boxpoint.vertex_min[0] <= (*root)->node_range_x[0] && boxpoint.vertex_max[0] > (*root)->node_range_x[1] && boxpoint.vertex_min[1] <= (*root)->node_range_y[0] && boxpoint.vertex_max[1] > (*root)->node_range_y[1] && boxpoint.vertex_min[2] <= (*root)->node_range_z[0] && boxpoint.vertex_max[2] > (*root)->node_range_z[1])
    {
        (*root)->tree_deleted = true;
        (*root)->point_deleted = true;
        (*root)->need_push_down_to_left = true;
        (*root)->need_push_down_to_right = true;
        tmp_counter = (*root)->tree_size - (*root)->invalid_point_num;
        (*root)->invalid_point_num = (*root)->tree_size;
        if (is_downsample)
        {
            (*root)->tree_downsample_deleted = true;
            (*root)->point_downsample_deleted = true;
            (*root)->down_del_num = (*root)->tree_size;
        }
        return tmp_counter;
    }
    if (!(*root)->point_deleted && boxpoint.vertex_min[0] <= (*root)->point.x && boxpoint.vertex_max[0] > (*root)->point.x && boxpoint.vertex_min[1] <= (*root)->point.y && boxpoint.vertex_max[1] > (*root)->point.y && boxpoint.vertex_min[2] <= (*root)->point.z && boxpoint.vertex_max[2] > (*root)->point.z)
    {
        (*root)->point_deleted = true;
        tmp_counter += 1;
        if (is_downsample)
            (*root)->point_downsample_deleted = true;
    }

    tmp_counter += DeleteByRange(&((*root)->left_son_ptr), boxpoint, allow_rebuild, is_downsample);
    tmp_counter += DeleteByRange(&((*root)->right_son_ptr), boxpoint, allow_rebuild, is_downsample);

    Update(*root);
    bool need_rebuild = allow_rebuild & CriterionCheck((*root));
    if (need_rebuild)
        Rebuild(root);
    if ((*root) != nullptr)
        (*root)->working_flag = false;
    return tmp_counter;
}

template <typename PointType>
void KD_TREE<PointType>::DeletePoints(PointVector& PointToDel) {
  for (int i = 0; i < PointToDel.size(); i++) {
    DeleteByPoint(&root_node_, PointToDel[i], true);
  }
  return;
}

template <typename PointType>
void KD_TREE<PointType>::AcquireRemovedPoints(PointVector& removed_points) {
  for (int i = 0; i < points_deleted_.size(); i++) {
    removed_points.push_back(points_deleted_[i]);
  }
  points_deleted_.clear();
  return;
}

template <typename PointType>
void KD_TREE<PointType>::BuildTree(KD_TREE_NODE** root, int l, int r,
                                   PointVector& Storage) {
  if (l > r) return;
  *root = new KD_TREE_NODE;
  InitTreeNode(*root);
  int mid = (l + r) >> 1;
  int div_axis = 0;
  int i;
  // Find the best division Axis
  float min_value[3] = {INFINITY, INFINITY, INFINITY};
  float max_value[3] = {-INFINITY, -INFINITY, -INFINITY};
  float dim_range[3] = {0, 0, 0};
  for (i = l; i <= r; i++) {
    min_value[0] = std::min(min_value[0], Storage[i].x);
    min_value[1] = std::min(min_value[1], Storage[i].y);
    min_value[2] = std::min(min_value[2], Storage[i].z);
    max_value[0] = std::max(max_value[0], Storage[i].x);
    max_value[1] = std::max(max_value[1], Storage[i].y);
    max_value[2] = std::max(max_value[2], Storage[i].z);
  }
  // Select the longest dimension as division axis
  for (i = 0; i < 3; i++) dim_range[i] = max_value[i] - min_value[i];
  for (i = 1; i < 3; i++)
    if (dim_range[i] > dim_range[div_axis]) div_axis = i;
  // Divide by the division axis and recursively build.

  (*root)->division_axis = div_axis;
  switch (div_axis) {
    case 0:
      std::nth_element(begin(Storage) + l, begin(Storage) + mid,
                  begin(Storage) + r + 1, PointCmpX);
      break;
    case 1:
      std::nth_element(begin(Storage) + l, begin(Storage) + mid,
                  begin(Storage) + r + 1, PointCmpY);
      break;
    case 2:
      std::nth_element(begin(Storage) + l, begin(Storage) + mid,
                  begin(Storage) + r + 1, PointCmpZ);
      break;
    default:
      std::nth_element(begin(Storage) + l, begin(Storage) + mid,
                  begin(Storage) + r + 1, PointCmpX);
      break;
  }
  (*root)->point = Storage[mid];
  KD_TREE_NODE *left_son = nullptr, *right_son = nullptr;
  BuildTree(&left_son, l, mid - 1, Storage);
  BuildTree(&right_son, mid + 1, r, Storage);
  (*root)->left_son_ptr = left_son;
  (*root)->right_son_ptr = right_son;
  Update((*root));
  return;
}

template <typename PointType>
void KD_TREE<PointType>::Rebuild(KD_TREE_NODE** root) {
  KD_TREE_NODE* father_ptr;
  father_ptr = (*root)->father_ptr;
  int size_rec = (*root)->tree_size;
  pcl_storage_.clear();
  Flatten(*root, pcl_storage_, DELETE_POINTS_REC);
  DeleteTreeNodes(root);
  BuildTree(root, 0, pcl_storage_.size() - 1, pcl_storage_);
  if (*root != nullptr) (*root)->father_ptr = father_ptr;
  if (*root == root_node_) static_root_node_->left_son_ptr = *root;
  return;
}

template <typename PointType>
void KD_TREE<PointType>::DeleteByPoint(KD_TREE_NODE** root, PointType point,
                                         bool allow_rebuild) {
  if ((*root) == nullptr || (*root)->tree_deleted) return;
  (*root)->working_flag = true;
  PushDown(*root);
  if (SamePoint((*root)->point, point) && !(*root)->point_deleted) {
    (*root)->point_deleted = true;
    (*root)->invalid_point_num += 1;
    if ((*root)->invalid_point_num == (*root)->tree_size)
      (*root)->tree_deleted = true;
    return;
  }

  if (((*root)->division_axis == 0 && point.x < (*root)->point.x) ||
      ((*root)->division_axis == 1 && point.y < (*root)->point.y) ||
      ((*root)->division_axis == 2 && point.z < (*root)->point.z)) {
    DeleteByPoint(&(*root)->left_son_ptr, point, allow_rebuild);
  } else {
    DeleteByPoint(&(*root)->right_son_ptr, point, allow_rebuild);
  }
  Update(*root);
  bool need_rebuild = allow_rebuild & CriterionCheck((*root));
  if (need_rebuild) Rebuild(root);
  if ((*root) != nullptr) (*root)->working_flag = false;
  return;
}

template <typename PointType>
void KD_TREE<PointType>::AddByPoint(KD_TREE_NODE** root, PointType point,
                                      bool allow_rebuild, int father_axis) {
  if (*root == nullptr) {
    *root = new KD_TREE_NODE;
    InitTreeNode(*root);
    (*root)->point = point;
    (*root)->division_axis = (father_axis + 1) % 3;
    Update(*root);
    return;
  }
  (*root)->working_flag = true;
  PushDown(*root);
  if (((*root)->division_axis == 0 && point.x < (*root)->point.x) ||
      ((*root)->division_axis == 1 && point.y < (*root)->point.y) ||
      ((*root)->division_axis == 2 && point.z < (*root)->point.z)) {
    AddByPoint(&(*root)->left_son_ptr, point, allow_rebuild,
                 (*root)->division_axis);
  } else {
    AddByPoint(&(*root)->right_son_ptr, point, allow_rebuild,
                 (*root)->division_axis);
  }
  Update(*root);
  bool need_rebuild = allow_rebuild & CriterionCheck((*root));
  if (need_rebuild) Rebuild(root);
  if ((*root) != nullptr) (*root)->working_flag = false;
  return;
}

template <typename PointType>
void KD_TREE<PointType>::Search(KD_TREE_NODE* root, int k_nearest,
                                PointType point, MANUAL_HEAP& q,
                                double max_dist) {
  if (root == nullptr || root->tree_deleted) return;
  double cur_dist = CalcBoxDist(root, point);
  double max_dist_sqr = max_dist * max_dist;
  if (cur_dist > max_dist_sqr) return;
  int retval;
  if (root->need_push_down_to_left || root->need_push_down_to_right) {
    retval = pthread_mutex_trylock(&(root->push_down_mutex_lock));
    if (retval == 0) {
      PushDown(root);
      pthread_mutex_unlock(&(root->push_down_mutex_lock));
    } else {
      pthread_mutex_lock(&(root->push_down_mutex_lock));
      pthread_mutex_unlock(&(root->push_down_mutex_lock));
    }
  }
  if (!root->point_deleted) {
    float dist = CalcDist(point, root->point);
    if (dist <= max_dist_sqr && (q.size() < k_nearest || dist < q.top().dist)) {
      if (q.size() >= k_nearest) q.pop();
      PointType_CMP current_point{root->point, dist};
      q.push(current_point);
    }
  }
  int cur_search_counter;
  float dist_left_node = CalcBoxDist(root->left_son_ptr, point);
  float dist_right_node = CalcBoxDist(root->right_son_ptr, point);
  if (q.size() < k_nearest ||
      dist_left_node < q.top().dist && dist_right_node < q.top().dist) {
    if (dist_left_node <= dist_right_node) {
      Search(root->left_son_ptr, k_nearest, point, q, max_dist);
      if (q.size() < k_nearest || dist_right_node < q.top().dist) {
        Search(root->right_son_ptr, k_nearest, point, q, max_dist);
      }
    } else {
      Search(root->right_son_ptr, k_nearest, point, q, max_dist);
      if (q.size() < k_nearest || dist_left_node < q.top().dist) {
        Search(root->left_son_ptr, k_nearest, point, q, max_dist);
      }
    }
  } else {
    if (dist_left_node < q.top().dist) {
      Search(root->left_son_ptr, k_nearest, point, q, max_dist);
    }
    if (dist_right_node < q.top().dist) {
      Search(root->right_son_ptr, k_nearest, point, q, max_dist);
    }
  }
  return;
}

template <typename PointType>
void KD_TREE<PointType>::SearchByRange(KD_TREE_NODE* root,
                                         BoxPointType boxpoint,
                                         PointVector& Storage) {
  if (root == nullptr) return;
  PushDown(root);
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
    Flatten(root, Storage, NOT_RECORD);
    return;
  }
  if (boxpoint.vertex_min[0] <= root->point.x &&
      boxpoint.vertex_max[0] > root->point.x &&
      boxpoint.vertex_min[1] <= root->point.y &&
      boxpoint.vertex_max[1] > root->point.y &&
      boxpoint.vertex_min[2] <= root->point.z &&
      boxpoint.vertex_max[2] > root->point.z) {
    if (!root->point_deleted) Storage.push_back(root->point);
  }
  SearchByRange(root->left_son_ptr, boxpoint, Storage);
  SearchByRange(root->right_son_ptr, boxpoint, Storage);
  return;
}

template <typename PointType>
void KD_TREE<PointType>::SearchByRadius(KD_TREE_NODE* root, PointType point,
                                          float radius, PointVector& Storage) {
  if (root == nullptr) return;
  PushDown(root);
  PointType range_center;
  range_center.x = (root->node_range_x[0] + root->node_range_x[1]) * 0.5;
  range_center.y = (root->node_range_y[0] + root->node_range_y[1]) * 0.5;
  range_center.z = (root->node_range_z[0] + root->node_range_z[1]) * 0.5;
  float dist = sqrt(CalcDist(range_center, point));
  if (dist > radius + sqrt(root->radius_sq)) return;
  if (dist <= radius - sqrt(root->radius_sq)) {
    Flatten(root, Storage, NOT_RECORD);
    return;
  }
  if (!root->point_deleted &&
      CalcDist(root->point, point) <= radius * radius) {
    Storage.push_back(root->point);
  }
  SearchByRadius(root->left_son_ptr, point, radius, Storage);
  SearchByRadius(root->right_son_ptr, point, radius, Storage);
  return;
}

template <typename PointType>
bool KD_TREE<PointType>::CriterionCheck(KD_TREE_NODE* root) {
  if (root->tree_size <= Minimal_Unbalanced_Tree_Size) {
    return false;
  }
  float balance_evaluation = 0.0f;
  float delete_evaluation = 0.0f;
  KD_TREE_NODE* son_ptr = root->left_son_ptr;
  if (son_ptr == nullptr) son_ptr = root->right_son_ptr;
  delete_evaluation = float(root->invalid_point_num) / root->tree_size;
  balance_evaluation = float(son_ptr->tree_size) / (root->tree_size - 1);
  if (delete_evaluation > delete_criterion_param_) {
    return true;
  }
  if (balance_evaluation > balance_criterion_param_ ||
      balance_evaluation < 1 - balance_criterion_param_) {
    return true;
  }
  return false;
}

template <typename PointType>
void KD_TREE<PointType>::PushDown(KD_TREE_NODE* root) {
  if (root == nullptr) return;
  if (root->need_push_down_to_left && root->left_son_ptr != nullptr) {
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
      root->left_son_ptr->down_del_num = root->left_son_ptr->tree_size;
    if (root->tree_deleted)
      root->left_son_ptr->invalid_point_num = root->left_son_ptr->tree_size;
    else
      root->left_son_ptr->invalid_point_num = root->left_son_ptr->down_del_num;
    root->left_son_ptr->need_push_down_to_left = true;
    root->left_son_ptr->need_push_down_to_right = true;
    root->need_push_down_to_left = false;
  }
  if (root->need_push_down_to_right && root->right_son_ptr != nullptr) {
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
      root->right_son_ptr->down_del_num = root->right_son_ptr->tree_size;
    if (root->tree_deleted)
      root->right_son_ptr->invalid_point_num = root->right_son_ptr->tree_size;
    else
      root->right_son_ptr->invalid_point_num =
          root->right_son_ptr->down_del_num;
    root->right_son_ptr->need_push_down_to_left = true;
    root->right_son_ptr->need_push_down_to_right = true;
    root->need_push_down_to_right = false;
  }
  return;
}

template <typename PointType>
void KD_TREE<PointType>::Update(KD_TREE_NODE* root) {
  KD_TREE_NODE* left_son_ptr = root->left_son_ptr;
  KD_TREE_NODE* right_son_ptr = root->right_son_ptr;
  float tmp_range_x[2] = {INFINITY, -INFINITY};
  float tmp_range_y[2] = {INFINITY, -INFINITY};
  float tmp_range_z[2] = {INFINITY, -INFINITY};
  // Update Tree Size
  if (left_son_ptr != nullptr && right_son_ptr != nullptr) {
    root->tree_size = left_son_ptr->tree_size + right_son_ptr->tree_size + 1;
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
         !root->point_deleted)) {
      tmp_range_x[0] = std::min(
          std::min(left_son_ptr->node_range_x[0], right_son_ptr->node_range_x[0]),
          root->point.x);
      tmp_range_x[1] = std::max(
          std::max(left_son_ptr->node_range_x[1], right_son_ptr->node_range_x[1]),
          root->point.x);
      tmp_range_y[0] = std::min(
          std::min(left_son_ptr->node_range_y[0], right_son_ptr->node_range_y[0]),
          root->point.y);
      tmp_range_y[1] = std::max(
          std::max(left_son_ptr->node_range_y[1], right_son_ptr->node_range_y[1]),
          root->point.y);
      tmp_range_z[0] = std::min(
          std::min(left_son_ptr->node_range_z[0], right_son_ptr->node_range_z[0]),
          root->point.z);
      tmp_range_z[1] = std::max(
          std::max(left_son_ptr->node_range_z[1], right_son_ptr->node_range_z[1]),
          root->point.z);
    } else {
      if (!left_son_ptr->tree_deleted) {
        tmp_range_x[0] = std::min(tmp_range_x[0], left_son_ptr->node_range_x[0]);
        tmp_range_x[1] = std::max(tmp_range_x[1], left_son_ptr->node_range_x[1]);
        tmp_range_y[0] = std::min(tmp_range_y[0], left_son_ptr->node_range_y[0]);
        tmp_range_y[1] = std::max(tmp_range_y[1], left_son_ptr->node_range_y[1]);
        tmp_range_z[0] = std::min(tmp_range_z[0], left_son_ptr->node_range_z[0]);
        tmp_range_z[1] = std::max(tmp_range_z[1], left_son_ptr->node_range_z[1]);
      }
      if (!right_son_ptr->tree_deleted) {
        tmp_range_x[0] = std::min(tmp_range_x[0], right_son_ptr->node_range_x[0]);
        tmp_range_x[1] = std::max(tmp_range_x[1], right_son_ptr->node_range_x[1]);
        tmp_range_y[0] = std::min(tmp_range_y[0], right_son_ptr->node_range_y[0]);
        tmp_range_y[1] = std::max(tmp_range_y[1], right_son_ptr->node_range_y[1]);
        tmp_range_z[0] = std::min(tmp_range_z[0], right_son_ptr->node_range_z[0]);
        tmp_range_z[1] = std::max(tmp_range_z[1], right_son_ptr->node_range_z[1]);
      }
      if (!root->point_deleted) {
        tmp_range_x[0] = std::min(tmp_range_x[0], root->point.x);
        tmp_range_x[1] = std::max(tmp_range_x[1], root->point.x);
        tmp_range_y[0] = std::min(tmp_range_y[0], root->point.y);
        tmp_range_y[1] = std::max(tmp_range_y[1], root->point.y);
        tmp_range_z[0] = std::min(tmp_range_z[0], root->point.z);
        tmp_range_z[1] = std::max(tmp_range_z[1], root->point.z);
      }
    }
  } else if (left_son_ptr != nullptr) {
    root->tree_size = left_son_ptr->tree_size + 1;
    root->invalid_point_num =
        left_son_ptr->invalid_point_num + (root->point_deleted ? 1 : 0);
    root->down_del_num =
        left_son_ptr->down_del_num + (root->point_downsample_deleted ? 1 : 0);
    root->tree_downsample_deleted =
        left_son_ptr->tree_downsample_deleted & root->point_downsample_deleted;
    root->tree_deleted = left_son_ptr->tree_deleted && root->point_deleted;
    if (root->tree_deleted ||
        (!left_son_ptr->tree_deleted && !root->point_deleted)) {
      tmp_range_x[0] = std::min(left_son_ptr->node_range_x[0], root->point.x);
      tmp_range_x[1] = std::max(left_son_ptr->node_range_x[1], root->point.x);
      tmp_range_y[0] = std::min(left_son_ptr->node_range_y[0], root->point.y);
      tmp_range_y[1] = std::max(left_son_ptr->node_range_y[1], root->point.y);
      tmp_range_z[0] = std::min(left_son_ptr->node_range_z[0], root->point.z);
      tmp_range_z[1] = std::max(left_son_ptr->node_range_z[1], root->point.z);
    } else {
      if (!left_son_ptr->tree_deleted) {
        tmp_range_x[0] = std::min(tmp_range_x[0], left_son_ptr->node_range_x[0]);
        tmp_range_x[1] = std::max(tmp_range_x[1], left_son_ptr->node_range_x[1]);
        tmp_range_y[0] = std::min(tmp_range_y[0], left_son_ptr->node_range_y[0]);
        tmp_range_y[1] = std::max(tmp_range_y[1], left_son_ptr->node_range_y[1]);
        tmp_range_z[0] = std::min(tmp_range_z[0], left_son_ptr->node_range_z[0]);
        tmp_range_z[1] = std::max(tmp_range_z[1], left_son_ptr->node_range_z[1]);
      }
      if (!root->point_deleted) {
        tmp_range_x[0] = std::min(tmp_range_x[0], root->point.x);
        tmp_range_x[1] = std::max(tmp_range_x[1], root->point.x);
        tmp_range_y[0] = std::min(tmp_range_y[0], root->point.y);
        tmp_range_y[1] = std::max(tmp_range_y[1], root->point.y);
        tmp_range_z[0] = std::min(tmp_range_z[0], root->point.z);
        tmp_range_z[1] = std::max(tmp_range_z[1], root->point.z);
      }
    }

  } else if (right_son_ptr != nullptr) {
    root->tree_size = right_son_ptr->tree_size + 1;
    root->invalid_point_num =
        right_son_ptr->invalid_point_num + (root->point_deleted ? 1 : 0);
    root->down_del_num =
        right_son_ptr->down_del_num + (root->point_downsample_deleted ? 1 : 0);
    root->tree_downsample_deleted =
        right_son_ptr->tree_downsample_deleted & root->point_downsample_deleted;
    root->tree_deleted = right_son_ptr->tree_deleted && root->point_deleted;
    if (root->tree_deleted ||
        (!right_son_ptr->tree_deleted && !root->point_deleted)) {
      tmp_range_x[0] = std::min(right_son_ptr->node_range_x[0], root->point.x);
      tmp_range_x[1] = std::max(right_son_ptr->node_range_x[1], root->point.x);
      tmp_range_y[0] = std::min(right_son_ptr->node_range_y[0], root->point.y);
      tmp_range_y[1] = std::max(right_son_ptr->node_range_y[1], root->point.y);
      tmp_range_z[0] = std::min(right_son_ptr->node_range_z[0], root->point.z);
      tmp_range_z[1] = std::max(right_son_ptr->node_range_z[1], root->point.z);
    } else {
      if (!right_son_ptr->tree_deleted) {
        tmp_range_x[0] = std::min(tmp_range_x[0], right_son_ptr->node_range_x[0]);
        tmp_range_x[1] = std::max(tmp_range_x[1], right_son_ptr->node_range_x[1]);
        tmp_range_y[0] = std::min(tmp_range_y[0], right_son_ptr->node_range_y[0]);
        tmp_range_y[1] = std::max(tmp_range_y[1], right_son_ptr->node_range_y[1]);
        tmp_range_z[0] = std::min(tmp_range_z[0], right_son_ptr->node_range_z[0]);
        tmp_range_z[1] = std::max(tmp_range_z[1], right_son_ptr->node_range_z[1]);
      }
      if (!root->point_deleted) {
        tmp_range_x[0] = std::min(tmp_range_x[0], root->point.x);
        tmp_range_x[1] = std::max(tmp_range_x[1], root->point.x);
        tmp_range_y[0] = std::min(tmp_range_y[0], root->point.y);
        tmp_range_y[1] = std::max(tmp_range_y[1], root->point.y);
        tmp_range_z[0] = std::min(tmp_range_z[0], root->point.z);
        tmp_range_z[1] = std::max(tmp_range_z[1], root->point.z);
      }
    }
  } else {
    root->tree_size = 1;
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
  root->radius_sq = x_L * x_L + y_L * y_L + z_L * z_L;
  if (left_son_ptr != nullptr) left_son_ptr->father_ptr = root;
  if (right_son_ptr != nullptr) right_son_ptr->father_ptr = root;
  if (root == root_node_ && root->tree_size > 3) {
    KD_TREE_NODE* son_ptr = root->left_son_ptr;
    if (son_ptr == nullptr) son_ptr = root->right_son_ptr;
    float tmp_bal = float(son_ptr->tree_size) / (root->tree_size - 1);
    root->alpha_del = float(root->invalid_point_num) / root->tree_size;
    root->alpha_bal = (tmp_bal >= 0.5 - EPSS) ? tmp_bal : 1 - tmp_bal;
  }
  return;
}

template <typename PointType>
void KD_TREE<PointType>::Flatten(KD_TREE_NODE *root, PointVector &Storage, DeletePointStorageSet storage_type)
{
    if (root == nullptr)
        return;
    PushDown(root);
    if (!root->point_deleted)
    {
        Storage.push_back(root->point);
    }
    Flatten(root->left_son_ptr, Storage, storage_type);
    Flatten(root->right_son_ptr, Storage, storage_type);
    switch (storage_type)
    {
    case NOT_RECORD:
        break;
    case DELETE_POINTS_REC:
        if (root->point_deleted && !root->point_downsample_deleted)
        {
            points_deleted_.push_back(root->point);
        }
        break;
    default:
        break;
    }
    return;
}

template <typename PointType>
void KD_TREE<PointType>::DeleteTreeNodes(KD_TREE_NODE** root) {
  if (*root == nullptr) return;
  PushDown(*root);
  DeleteTreeNodes(&(*root)->left_son_ptr);
  DeleteTreeNodes(&(*root)->right_son_ptr);

  pthread_mutex_destroy(&(*root)->push_down_mutex_lock);
  delete *root;
  *root = nullptr;

  return;
}

template <typename PointType>
bool KD_TREE<PointType>::SamePoint(PointType a, PointType b) {
  return (fabs(a.x - b.x) < EPSS && fabs(a.y - b.y) < EPSS &&
          fabs(a.z - b.z) < EPSS);
}

template <typename PointType>
float KD_TREE<PointType>::CalcDist(PointType a, PointType b) {
  float dist = 0.0f;
  dist = (a.x - b.x) * (a.x - b.x) + (a.y - b.y) * (a.y - b.y) +
         (a.z - b.z) * (a.z - b.z);
  return dist;
}

template <typename PointType>
float KD_TREE<PointType>::CalcBoxDist(KD_TREE_NODE* node, PointType point) {
  if (node == nullptr) return INFINITY;
  float min_dist = 0.0;
  if (point.x < node->node_range_x[0])
    min_dist +=
        (point.x - node->node_range_x[0]) * (point.x - node->node_range_x[0]);
  if (point.x > node->node_range_x[1])
    min_dist +=
        (point.x - node->node_range_x[1]) * (point.x - node->node_range_x[1]);
  if (point.y < node->node_range_y[0])
    min_dist +=
        (point.y - node->node_range_y[0]) * (point.y - node->node_range_y[0]);
  if (point.y > node->node_range_y[1])
    min_dist +=
        (point.y - node->node_range_y[1]) * (point.y - node->node_range_y[1]);
  if (point.z < node->node_range_z[0])
    min_dist +=
        (point.z - node->node_range_z[0]) * (point.z - node->node_range_z[0]);
  if (point.z > node->node_range_z[1])
    min_dist +=
        (point.z - node->node_range_z[1]) * (point.z - node->node_range_z[1]);
  return min_dist;
}

template <typename PointType>
bool KD_TREE<PointType>::PointCmpX(PointType a, PointType b) {
  return a.x < b.x;
}
template <typename PointType>
bool KD_TREE<PointType>::PointCmpY(PointType a, PointType b) {
  return a.y < b.y;
}
template <typename PointType>
bool KD_TREE<PointType>::PointCmpZ(PointType a, PointType b) {
  return a.z < b.z;
}

// manual queue
template <typename T>
void MANUAL_Q<T>::clear() {
  head = 0;
  tail = 0;
  counter = 0;
  is_empty = true;
  return;
}

template <typename T>
void MANUAL_Q<T>::pop() {
  if (counter == 0) return;
  head++;
  head %= Q_LEN;
  counter--;
  if (counter == 0) is_empty = true;
  return;
}

template <typename T>
T MANUAL_Q<T>::front() {
  return q[head];
}

template <typename T>
T MANUAL_Q<T>::back() {
  return q[tail];
}

template <typename T>
void MANUAL_Q<T>::push(T op) {
  q[tail] = op;
  counter++;
  if (is_empty) is_empty = false;
  tail++;
  tail %= Q_LEN;
}

template <typename T>
bool MANUAL_Q<T>::empty() {
  return is_empty;
}

template <typename T>
int MANUAL_Q<T>::size() {
  return counter;
}

template class KD_TREE<pcl::PointXYZ>;
template class KD_TREE<pcl::PointXYZI>;
template class KD_TREE<pcl::PointXYZINormal>;