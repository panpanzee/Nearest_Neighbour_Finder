#include "kd_tree.h"
#include <boost/lexical_cast.hpp>
#include <boost/bind.hpp>
#include <boost/function.hpp>
#include <iostream>
#include <ostream>
#include <queue>
#include <math.h>
#include <algorithm>
#include <limits>
#include <stdexcept>
#include <boost/date_time.hpp>
#include <boost/thread.hpp> 
#include <numeric>

/****************
Core implemention of kd_tree.
I adopt the standard design of kd_tree here with several implmentatino tricks to boost the performance. Because the progarm is consider to be running on a relatively large data set, I care both the space complesity and the time complexity.
1. Using a vector of shared pointer to keep track of the points, such that the actual is stored only once.
2. Using iterators of shared pointers for sorting and partitioning, such that the actual data is never copied around.
3. Checking the actual running time in searching, and early stop if 1 sec limit is reached, to make sure we always finish within 1 sec.
4. Supporting quick_select algorithm for finding median and partitionion around it. It is O(N) in average, compared with O(nlogn) for sorting and find median.
5. Normalising the data before constructing the kd_tree and providing support for weighted euclidean distance.
6. Leaf node size is customizable, to allow the client code to optimise the performance.
7.Supporting fast but not optimised tree constrution algorithm. This is usefull when the data set is too big and the launching time takes so long.
***************/


using namespace std;
using namespace boost;

typedef vector<shared_ptr<KdPointConvertable> >::iterator KDPointIterator;
typedef vector<shared_ptr<KdPointConvertable> >::const_iterator KDPointConstIterator;

namespace {

  size_t get_best_partition_dimension(const KDPointConstIterator begin,
                                      const KDPointConstIterator end)
  {
    const size_t num_of_dimensions = (*begin)->get_dimensions();
    vector<double> sums(num_of_dimensions, 0);
    vector<double> variances(num_of_dimensions, 0);
    for(KDPointConstIterator it = begin; 
        it != end; it++) {
      if((*it)->get_dimensions() != num_of_dimensions)
        throw runtime_error("Dimensionality not the same.");

      for(int dimension = 0; dimension < num_of_dimensions; dimension++) {
        double value = (*it)->get_cordinates(dimension);
        sums[dimension] += value;
      }
    }

    size_t data_set_size = end - begin;
    for(KDPointConstIterator it = begin; 
        it != end; it++) {
      for(int dimension = 0; dimension < num_of_dimensions; dimension++) {
        double value = (*it)->get_cordinates(dimension);
        variances[dimension] += pow(value - sums[dimension]/data_set_size, 2);
      }
    }

    return distance(variances.begin(), max_element(variances.begin(), variances.end()));
  }

  size_t get_next_partition_dimension(const KDPointConstIterator begin,
                                      const KDPointConstIterator end)
  {
    static size_t partition_dimension = -1;
    partition_dimension++;
    const size_t num_of_dimensions = (*begin)->get_dimensions();
    if(partition_dimension == num_of_dimensions)
      partition_dimension = 0;
    return partition_dimension;
  }

  bool compare(const shared_ptr<KdPointConvertable>& lhs, const shared_ptr<KdPointConvertable>& rhs, const size_t dimension)
  {
    return lhs->get_cordinates(dimension) <= rhs->get_cordinates(dimension);
  }

  KDPointIterator partition(const KDPointIterator begin, const KDPointIterator end, const KDPointIterator pivot_iterator, const size_t partition_dimension)
  {
    shared_ptr<KdPointConvertable> pivot_value = *pivot_iterator;
    iter_swap(pivot_iterator, end-1);
    KDPointIterator split_iterarot = begin;
    for(KDPointIterator it = begin; it != end - 2; it++){
      if(compare(*it, pivot_value, partition_dimension)){
        iter_swap(split_iterarot, it);
        split_iterarot++;
      }
    }
    iter_swap(split_iterarot, end-1);
    return split_iterarot;
  }

  double get_any_value_and_partition(const KDPointIterator in_begin, 
                                     const KDPointIterator in_end,
                                     const size_t in_partition_dimension,
                                     KDPointIterator& out_begin_of_first_partition)
  {
    const KDPointIterator split_iterator = partition(in_begin, in_end, in_begin, in_partition_dimension);
    out_begin_of_first_partition = split_iterator + 1;
    return (*split_iterator)->get_cordinates(in_partition_dimension);
  }

  double get_median_and_partition_by_quickselect(const KDPointIterator in_begin, 
                                                 const KDPointIterator in_end,
                                                 const size_t in_partition_dimension,
                                                 KDPointIterator& out_begin_of_first_partition)
  {
    size_t length = in_end - in_begin;
    size_t target_index = length%2 == 0 ? length/2 - 1: length/2;
    KDPointIterator begin = in_begin;
    KDPointIterator end = in_end;
    while(begin < end){
      KDPointIterator split_iterator = partition(begin, end, begin, in_partition_dimension);
      const size_t split_point_index = split_iterator - begin;
      if(split_point_index == target_index){
        out_begin_of_first_partition = split_iterator + 1;
        return (*split_iterator)->get_cordinates(in_partition_dimension);
      }

      if(split_point_index < target_index){
        begin = split_iterator + 1;
        target_index = target_index - split_point_index - 1;
      }else{
        end = split_iterator + 1;
      }     
    }
  }

  double get_median_and_partition_by_sorting(const KDPointIterator in_begin, 
                                             const KDPointIterator in_end,
                                             const size_t in_partition_dimension,
                                             KDPointIterator& out_begin_of_first_partition)
  {
    sort(in_begin, in_end, boost::bind(compare, _1, _2, in_partition_dimension));
    const size_t length = in_end - in_begin;
    out_begin_of_first_partition = in_begin + length/2;
    double median_value = (*out_begin_of_first_partition)->get_cordinates(in_partition_dimension);
    double rval;
    if(length % 2 == 0){
      rval = (median_value + (*(out_begin_of_first_partition-1))->get_cordinates(in_partition_dimension)) / 2;
    }else{
      rval = median_value;
    }

    KDPointIterator next = out_begin_of_first_partition + 1;
    while((next != in_end) && (median_value == (*next)->get_cordinates(in_partition_dimension))){
      out_begin_of_first_partition++;
      next = out_begin_of_first_partition + 1;
    }
    return rval;
  }
}

const int Normaliser::g_normalisation_range = 10000; 

double Normaliser::calculate_new_value(const double old_value, const double mean, const double range)
{
  if(fabs((old_value - mean) / range) >1 )
    int a = 0;
  return (old_value - mean) / range * g_normalisation_range;
}

void Normaliser::transform_to_normalisation_space(KdPointConvertable* point)
{
  const size_t num_of_dimensions = point->get_dimensions();
  for(int dimension = 0; dimension < num_of_dimensions; dimension++) {
    double old_value = point->get_cordinates(dimension);
    double new_value = calculate_new_value(old_value, m_means[dimension], m_ranges[dimension]);
    point->set_cordinates(dimension, new_value);
  }  
}

void Normaliser::normalisation(vector<shared_ptr<KdPointConvertable> >& points)
  {
    if(points.empty())
      return;

    const size_t num_of_dimensions = points[0]->get_dimensions();
    vector<double> sums(num_of_dimensions, 0);
    m_means.resize(num_of_dimensions, 0);
    m_ranges.resize(num_of_dimensions, 0);
    vector<double> maxes(num_of_dimensions, numeric_limits<double>::min());
    vector<double> mins(num_of_dimensions, numeric_limits<double>::max());

    for(KDPointConstIterator it = points.begin(); 
        it != points.end(); it++) {
      if((*it)->get_dimensions() != num_of_dimensions)
        throw runtime_error("Dimensionality not the same.");

      for(int dimension = 0; dimension < num_of_dimensions; dimension++) {
        double value = (*it)->get_cordinates(dimension);
        sums[dimension] += value;
        if(value > maxes[dimension])
          maxes[dimension] = value;
        if(value < mins[dimension])
          mins[dimension] = value;
      }
    }

    size_t data_set_size = points.size();

    for(KDPointIterator it = points.begin(); 
        it != points.end(); it++) {
      for(int dimension = 0; dimension < num_of_dimensions; dimension++) {
        double old_value = (*it)->get_cordinates(dimension);
        m_means[dimension] = sums[dimension]/data_set_size;
        m_ranges[dimension] = maxes[dimension] - mins[dimension];
        double new_value = calculate_new_value(old_value, m_means[dimension], m_ranges[dimension]);
        (*it)->set_cordinates(dimension, new_value);
      }
    }    
  }

double KdPointConvertable::get_cordinates(const size_t axis) const
{
  if(axis >= m_cordinates.size())
    throw runtime_error(string("Dimension " + boost::lexical_cast<string>(axis) + " out of bound.").c_str());

  return m_cordinates[axis];
}

void KdPointConvertable::set_cordinates(const size_t axis, const double value)
{
  if(axis >= m_cordinates.size())
    throw runtime_error(string("Dimension " + boost::lexical_cast<string>(axis) + " out of bound.").c_str());

  m_cordinates[axis] = value;
}

namespace Distance{
vector<double> g_weights;
bool is_using_weights = false;
}
void Distance::set_weights(std::vector<double> weights)
{
  if(weights.empty())
    return;

  is_using_weights = true;
  double sum = accumulate(weights.begin(),weights.end(),0);
  for(int i = 0; i < weights.size(); i++){
    if(weights[i] <= 0)
      throw runtime_error("Weights have to be positive numbers.");
    g_weights.push_back(weights[i] / sum);
  }
}

double Distance::get_euclidean_distance(const KdPointConvertable* p1, const KdPointConvertable* p2)
{
  size_t dimension_p1 = p1->get_dimensions();
  size_t dimension_p2 = p2->get_dimensions();
  if( dimension_p1 != dimension_p2)
    throw runtime_error("Dimensionality doesn't match.");

  if(is_using_weights){
    size_t dimension_weight = g_weights.size();  
    if(dimension_weight != dimension_p1){
      throw runtime_error("Dimensionality doesn't match.");
    }
  }

  double sum = 0;
  for(int i = 0; i < dimension_p1; i++){
    double distance_in_dimension_i = p1->get_cordinates(i) - p2->get_cordinates(i);
    if(is_using_weights)
      distance_in_dimension_i *= g_weights[i];
    sum += pow(distance_in_dimension_i, 2);
  }

  return sum;
}

struct Candidate {
  double m_distance;
  shared_ptr<KdPointConvertable> m_content;
  Candidate(const double ditance, const shared_ptr<KdPointConvertable> content)
    : m_distance(ditance), m_content(content)
  {}

  bool operator<(const Candidate& rhs) const
  {
    return m_distance < rhs.m_distance;
  }
};

typedef priority_queue<Candidate, vector<Candidate>,less<vector<Candidate>::value_type> > candidates_queue;


namespace{
  void try_to_add_a_candidate(candidates_queue& current_candidates,
                              const size_t max_candidates_number,
                              const shared_ptr<KdPointConvertable> raw_candidate,
                              const KdPointConvertable* query_point)
  {
    double distance = Distance::get_euclidean_distance(raw_candidate.get(), query_point);
    if(current_candidates.size() < max_candidates_number){
      current_candidates.push(Candidate(distance, raw_candidate));
      return;
    }

    const Candidate& current_worst_candidate = current_candidates.top();
    if(current_worst_candidate.m_distance > distance){
      current_candidates.pop();
      current_candidates.push(Candidate(distance, raw_candidate));
    }
  }

  void try_to_add_candidates(candidates_queue& current_candidates,
                            const size_t max_candidates_number,
                            const KDPointConstIterator raw_candidate_begin,
                            const KDPointConstIterator raw_candidate_end,
                            const KdPointConvertable* query_point)
  {
    for(KDPointConstIterator it = raw_candidate_begin; it != raw_candidate_end; it++)
      try_to_add_a_candidate(current_candidates, max_candidates_number, *it, query_point);
  }
}

class KdNode {
public:
  KdNode(const KDPointIterator begin,
         const KDPointIterator end);
  bool is_leaf_node() const { return m_split_dimension == -1; }
  void get_NN(const size_t number_of_NN, const KdPointConvertable* query_point, candidates_queue& current_candidates, boost::function0<bool> is_time_out) const;
  static size_t s_leaf_node_size;
  static boost::function<size_t(const KDPointConstIterator, const KDPointConstIterator)> s_get_partition_dimension;
  static boost::function<double(const KDPointIterator, const KDPointIterator, const size_t, KDPointIterator&)> s_get_median_and_partition;
private:
  double m_split_value;
  size_t m_split_dimension;
  shared_ptr<KdNode> m_left_subtree;
  shared_ptr<KdNode> m_right_subtree;
  KDPointConstIterator m_leaf_node_content_begin;
  KDPointConstIterator m_leaf_node_content_end;
};

size_t KdNode::s_leaf_node_size = 1;
boost::function<size_t(const KDPointConstIterator, const KDPointConstIterator)> KdNode::s_get_partition_dimension = get_best_partition_dimension;
boost::function<double(const KDPointIterator, const KDPointIterator, const size_t, KDPointIterator&)> KdNode::s_get_median_and_partition = get_median_and_partition_by_sorting;

KdNode::KdNode(const KDPointIterator begin, const KDPointIterator end)
  : m_split_value(0), m_split_dimension(0)
{
  if(end <= begin)
    throw runtime_error("Ending iterator comes before begining iterator. Something went wrong.");

  if(end - begin <= s_leaf_node_size){
    m_split_dimension = -1;
    m_leaf_node_content_begin = begin;
    m_leaf_node_content_end = end;
    return;
  }

  m_split_dimension = s_get_partition_dimension(begin, end);
  KDPointIterator end_of_first_partition;
  m_split_value = s_get_median_and_partition(begin, end, m_split_dimension, end_of_first_partition);
  m_left_subtree.reset(new KdNode(begin, end_of_first_partition));
  m_right_subtree.reset(new KdNode(end_of_first_partition, end));
}

void KdNode::get_NN(const size_t number_of_NN, const KdPointConvertable* query_point, candidates_queue& current_candidates, boost::function0<bool> is_time_out) const
{
  if(is_time_out())
    return;

  if(is_leaf_node()){
    try_to_add_candidates(current_candidates, number_of_NN, m_leaf_node_content_begin, m_leaf_node_content_end, query_point);
    return;
  }

  bool search_left_subtree_first;
  if(query_point->get_cordinates(m_split_dimension) <= m_split_value){
    search_left_subtree_first = true;
  }else{
    search_left_subtree_first = false;
  }

  if(search_left_subtree_first)
    m_left_subtree->get_NN(number_of_NN, query_point, current_candidates, is_time_out);
  else
    m_right_subtree->get_NN(number_of_NN, query_point, current_candidates, is_time_out);

  const Candidate& current_worst_candidate = current_candidates.top();
  double distance_to_split_line = fabs(query_point->get_cordinates(m_split_dimension) - m_split_value);

  if(Distance::is_using_weights)
    distance_to_split_line *= Distance::g_weights[m_split_dimension];

  if((current_candidates.size() < number_of_NN) || (distance_to_split_line < current_worst_candidate.m_distance)){
    if(search_left_subtree_first)
      m_right_subtree->get_NN(number_of_NN, query_point, current_candidates, is_time_out); 
    else
      m_left_subtree->get_NN(number_of_NN, query_point, current_candidates, is_time_out); 
  }
}


KdTree::KdTree(vector<shared_ptr<KdPointConvertable> >& points, boost::shared_ptr<Normaliser> normaliser)
  : m_normaliser(normaliser)
{ 
  m_normaliser->normalisation(points);
  m_root.reset(new KdNode(points.begin(), points.end()));
}

void KdTree::set_leaf_node_size(const size_t size)
{ 
  KdNode::s_leaf_node_size = size; 
}

void KdTree::set_get_partition_dimension_algorithm(const GetPartitionDimensionAlgorithm algo)
{
  switch(algo){
  case LARGEST_VARIANCE:
    KdNode::s_get_partition_dimension = get_best_partition_dimension;
    break;
  case ROUND_ROBIN:
    KdNode::s_get_partition_dimension = get_next_partition_dimension;
    break;
  default:
    throw runtime_error("get_partition_dimension_algorithm not supported.");
  }
}

void KdTree::set_partition_algorithm(const PartitionAlgorithm algo)
{
  switch(algo){
  case SORTING:
    KdNode::s_get_median_and_partition = get_median_and_partition_by_sorting;
    break;
  case RANDOM:
    KdNode::s_get_median_and_partition = get_any_value_and_partition;
    break;
  case QUICK_SELECT:
    KdNode::s_get_median_and_partition = get_median_and_partition_by_quickselect;
    break;
  default:
    throw runtime_error("get_partition_dimension_algorithm not supported.");
  }
}

namespace {
  bool check_time_out(const struct timeval& begin)
  {
    struct timeval end;
    gettimeofday(&end, NULL);
    static double const microSecPerSec = 1000000;
    const double sec = (end.tv_sec + end.tv_usec/microSecPerSec) - (begin.tv_sec+begin.tv_usec/microSecPerSec);
    return sec > 1;
  }
}

vector<shared_ptr<KdPointConvertable> > KdTree::get_NN(const size_t number_of_NN, const KdPointConvertable* query_point) const
{
  struct timeval begin;
  gettimeofday(&begin, NULL);
  KdPointConvertable query_point_copy(*query_point);
  m_normaliser->transform_to_normalisation_space(&query_point_copy);
  candidates_queue candidates;
  m_root->get_NN(number_of_NN, &query_point_copy, candidates, boost::bind(check_time_out, begin));
  vector<shared_ptr<KdPointConvertable> > final_candidates;
  final_candidates.reserve(number_of_NN);
  while(!candidates.empty()){
    final_candidates.push_back(candidates.top().m_content);
    candidates.pop();
  }
  reverse(final_candidates.begin(), final_candidates.end());
  return final_candidates;
}
