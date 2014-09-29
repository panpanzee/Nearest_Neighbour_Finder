#include "stdlib.h"
#include <boost/scoped_ptr.hpp>
#include <boost/shared_ptr.hpp>
#include <vector>
#include <string>

// Kd tree node base class
// We use this to store the cordinates
// Derived class should populate the cordinates vector
// The dimension is determined by the actual size of the cordinates vector
class KdPointConvertable {
public:
  virtual void print(std::ostream& out) const {}
  double get_cordinates(const size_t axis) const;
  void set_cordinates(const size_t axis, const double value);
  size_t get_dimensions() const { return m_cordinates.size(); }
  
  KdPointConvertable()
  {}

protected:
  KdPointConvertable(const std::vector<double>& cordinates)
    :m_cordinates(cordinates)
  {}
  
  std::vector<double> m_cordinates;
};

// Helper functions to calculate the distance
// and set weights on dimensions
namespace Distance{
  double get_euclidean_distance(const KdPointConvertable* p1, const KdPointConvertable* p2);
  void set_weights(std::vector<double> weights);
}


// Normalisation class
// Used to normalise input data
// and also transform query points to normalised space
class Normaliser {
public:
  void normalisation(std::vector<boost::shared_ptr<KdPointConvertable> >& points);
  void transform_to_normalisation_space(KdPointConvertable* point);
private:
  inline double calculate_new_value(const double old_value, const double mean, const double range);
  const static int g_normalisation_range;
  std::vector<double> m_means;
  std::vector<double> m_ranges;
};

// Forward declaration of the Kd tree node class
class KdNode;

// Kd tree class
// Construct it with a vector of pointers to KdPointConvertable, and a normaliser.
class KdTree {
public:
  KdTree(std::vector<boost::shared_ptr<KdPointConvertable> >& points, boost::shared_ptr<Normaliser> normaliser);
  std::vector<boost::shared_ptr<KdPointConvertable> > get_NN(const size_t number_of_NN, const KdPointConvertable* query) const;
  static void set_leaf_node_size(const size_t size);

  enum GetPartitionDimensionAlgorithm{
    ROUND_ROBIN = 0,
    LARGEST_VARIANCE
  };  
  static void set_get_partition_dimension_algorithm(const GetPartitionDimensionAlgorithm algo);

  enum PartitionAlgorithm{
    SORTING = 0,
    RANDOM,
    QUICK_SELECT
  };

  static void set_partition_algorithm(const PartitionAlgorithm algo);
private:
  boost::shared_ptr<KdNode> m_root;
  boost::shared_ptr<Normaliser> m_normaliser;
};
