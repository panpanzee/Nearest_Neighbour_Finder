#include "kd_tree.h"
#include <boost/lexical_cast.hpp>
#include <boost/bind.hpp>
#include <iostream>     // cout, endl
#include <stdexcept>
#include <sys/time.h>

/*******************
Testing prog

The result of Kd_tree searching is compared with the result of brute force searching.
Upon launching, it generates a set of random points given the command line input and
use the points to construct the kd_tree and the brute force searcher.

There are 2 sets of parameters that can be specified through the command line, problem specification parameters and performance tweaking parameters.
The former includes the dimensionality of the points, the size of the input data set, the size of the test data set, and the number of nearest neighbors it's aiming to find. 
The latter contains the size of the kd_tree leaf node and 2 algorithm involved in kd_tree construction.
One is the algorithm used in partition selection. It can be selected between ROUND_BORIN and MOST_SPREAD_OUT( LARGEST_VARIANCE).
The other one is the partition algorithm. It can be selected among SORTING_AND_CHOOSE_MEDIAN, RANDOM, AND QUICK_SELECT_FOR_MEDIAN.

Then, it compare the result of each test point between kd_tree and brute force searching.
I only consider a test case to pass if all of the N nearest neighbors returned from both algorithms are the same.
The final pass rate of running time in each section are printed to the standard output.

********************/


using namespace std;
using namespace boost;

namespace Test{
  class KdPointTest : public KdPointConvertable {
  public:
    KdPointTest(const std::vector<double>& cordinates)
      : KdPointConvertable(cordinates)
    {}
  };

  bool compare(const shared_ptr<KdPointConvertable>& p1, const shared_ptr<KdPointConvertable>& p2, const KdPointConvertable* target_point)
  {
    double distance_p1 = Distance::get_euclidean_distance(p1.get(), target_point);
    double distance_p2 = Distance::get_euclidean_distance(p2.get(), target_point);
    return distance_p1 < distance_p2;
  } 
}

class BruteForce {
public:
  BruteForce(std::vector<boost::shared_ptr<KdPointConvertable> >& points, boost::shared_ptr<Normaliser> normaliser)
    :m_points(points), m_normaliser(normaliser)
  {}

  std::vector<boost::shared_ptr<KdPointConvertable> > get_NN(const size_t number_of_NN, const KdPointConvertable* query_point)
  {
    KdPointConvertable query_point_copy(*query_point);
    m_normaliser->transform_to_normalisation_space(&query_point_copy);

    vector<shared_ptr<KdPointConvertable> > current_candidates;
    for(int i = 0; i < m_points.size(); i++){
      if(current_candidates.size() < number_of_NN){
        current_candidates.push_back(m_points[i]);
        sort(current_candidates.begin(), current_candidates.end(), boost::bind(Test::compare, _1, _2, &query_point_copy));
        continue;
      }
        
      const KdPointConvertable* worst_candidate = current_candidates[number_of_NN-1].get();
      if(Distance::get_euclidean_distance(worst_candidate, &query_point_copy) > Distance::get_euclidean_distance(m_points[i].get(), &query_point_copy)){
          current_candidates[number_of_NN-1] = m_points[i];
          sort(current_candidates.begin(), current_candidates.end(), boost::bind(Test::compare, _1, _2, &query_point_copy));
      }
    }

    return current_candidates;
  }

private:
  vector<shared_ptr<KdPointConvertable> > m_points;
  boost::shared_ptr<Normaliser> m_normaliser;
};

namespace Test{
  boost::shared_ptr<KdPointConvertable> generate_a_point(const size_t dimensions)
  {
    vector<double> cordinates;
    for(int j = 0; j < dimensions; j++){
      cordinates.push_back(static_cast<double>(rand()) / static_cast<double>(RAND_MAX) * 1000);
    }
    return boost::shared_ptr<KdPointConvertable>(new KdPointTest(cordinates));
  }

  vector<boost::shared_ptr<KdPointConvertable> > generate_points(const size_t dimensions, const size_t data_set_size){
    vector<boost::shared_ptr<KdPointConvertable> > rval;
    for(int i = 0; i < data_set_size; i++){
      rval.push_back(generate_a_point(dimensions));
    }
    return rval;
  }

  template <class T>
  bool assert_equal(const T a, const T b)
  {
    if(a!=b){
      string error_msg(boost::lexical_cast<string>(a) + " != " + boost::lexical_cast<string>(b));
      cout<<error_msg<<endl;
      return false;
      //throw runtime_error(error_msg.c_str());
    }
    return true;
  }

  void print_running_time(const struct timeval& begin, const struct timeval& end, const string& algorithm)
  {
  static double const microSecPerSec = 1000000;

  const double sec = (end.tv_sec + end.tv_usec/microSecPerSec) - (begin.tv_sec+begin.tv_usec/microSecPerSec);
  cout<<algorithm<<" running time: "<<sec<<" secondes"<<endl<<flush;
  }
}



int main(int argc, char **argv)
{
  
  if(argc < 5){
    cout<<"Usage: NN_test <dimensions> <data_set_size> <number_of_test_points> <number_of_NN> [leaf_node_size (default=1)] [choose_partition_dimension_algolrithm, 0: ROUND_ROBIN, 1: LARGEST_VARIANCE (default=1)] [partition_algorithm, 0: SORTING, 1: RANDOM, 2:QUICK_SELECT (default=2)]"<<endl;
    return 1;
  }

  const int dimensions = atoi(argv[1]);
  const int data_set_size = atoi(argv[2]);
  const int number_of_test_points = atoi(argv[3]);
  const int number_of_NN = atoi(argv[4]);
  const int leaf_node_size = argc > 5 ? atoi(argv[5]) : 1;
  const int get_partition_dimension_algo = argc > 6 ? atoi(argv[6]) : 0;
  const int partition_algo = argc > 7 ? atoi(argv[7]) : 2;

  std::vector<boost::shared_ptr<KdPointConvertable> > sample_points(Test::generate_points(dimensions, data_set_size));

  struct timeval begin, end;
  boost::shared_ptr<Normaliser> normaliser(new Normaliser());

  gettimeofday(&begin, NULL);
  KdTree::set_leaf_node_size(leaf_node_size);

  if(get_partition_dimension_algo == 1)
    KdTree::set_get_partition_dimension_algorithm(KdTree::ROUND_ROBIN);
  if(partition_algo == 1)
    KdTree::set_partition_algorithm(KdTree::RANDOM);
  else if(partition_algo == 2)
    KdTree::set_partition_algorithm(KdTree::QUICK_SELECT);

  KdTree kd_tree(sample_points, normaliser);
  gettimeofday(&end, NULL);
  Test::print_running_time(begin, end, "KDTree construction");
  
  gettimeofday(&begin, NULL);
  BruteForce brute_force(sample_points, normaliser);
  gettimeofday(&end, NULL);
  Test::print_running_time(begin, end, "Brute force construction");

  int number_of_failed_test = 0;
  for(int j =0; j < number_of_test_points; j++){
    boost::shared_ptr<KdPointConvertable> target_point(Test::generate_a_point(dimensions));
    gettimeofday(&begin, NULL);
    const vector<boost::shared_ptr<KdPointConvertable> > from_Kd_tree(kd_tree.get_NN(number_of_NN, target_point.get()));
    gettimeofday(&end, NULL);
    Test::print_running_time(begin, end, "KDTree searching");

    gettimeofday(&begin, NULL);
    const vector<boost::shared_ptr<KdPointConvertable> > from_brute_force(brute_force.get_NN(number_of_NN, target_point.get()));
    gettimeofday(&end, NULL);
    Test::print_running_time(begin, end, "Brute force searching");

    if(!Test::assert_equal<size_t>(from_Kd_tree.size(), from_brute_force.size())){
      cout<<"Test point " << j+1 << " failed"<<endl<<endl;
      number_of_failed_test ++;
      continue;
    }

    int pass = 0;
    for(int i = 0; i < from_Kd_tree.size(); i++){
      if(Test::assert_equal<double>(Distance::get_euclidean_distance(from_Kd_tree[i].get(), target_point.get()), Distance::get_euclidean_distance(from_brute_force[i].get(), target_point.get())))
        pass++;
    }
    if(pass == from_Kd_tree.size()){
      cout<<"Test point " << j+1 << " passed"<<endl<<endl;
    }
    else{
      number_of_failed_test ++;
      cout<<pass << " out of " <<from_Kd_tree.size() << "neighbours pass for test point " << j+1 <<endl<<endl;
    }
  }

  if(number_of_failed_test == 0)
    cout<<"All test passed."<<endl;
  else
    cout<<number_of_failed_test << " out of " << number_of_test_points << "test failed."<<endl;
}
