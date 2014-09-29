#include <iostream>     // cout, endl
#include <fstream>      // fstream
#include <stdexcept>
#include <string.h>
#include <boost/lexical_cast.hpp>
#include<boost/tokenizer.hpp>
#include "kd_tree.h"


using namespace std;
using namespace boost;

class Person : public KdPointConvertable{
public:
  Person(const int age, const double latitude, const double longitude, const std::string& name);
  virtual void print(ostream& out) const; 
private:
  int m_age;
  double m_latitude;
  double m_longitude;
  string m_name;
};

Person::Person(const int age, const double latitude, const double longitude, const string& name)
  : m_age(age), m_latitude(latitude), m_longitude(longitude), m_name(name)
{
  m_cordinates.push_back(age);
  m_cordinates.push_back(latitude);
  m_cordinates.push_back(longitude);
}

void Person::print(ostream& out) const
{
  out << m_name<<endl;
  out<< "Age: "
      << m_age
      << ", Latitude: "
      << m_latitude
      << ", Longtitude: "
      << m_longitude
     << "."<<endl<<endl;
}

class PeopleCreater {
public:
  PeopleCreater(const string& file_names = "names.csv");
  boost::shared_ptr<KdPointConvertable> generate_a_person() const;
  vector<boost::shared_ptr<KdPointConvertable> > generate_people(const size_t number_of_people) const;
private:
  vector<string> m_first_names;
  vector<string> m_family_names;
  int m_min_age;
  int m_max_age;
  double m_min_latitude;
  double m_max_latitude;
  double m_min_longtitude;
  double m_max_longtitude;
};

PeopleCreater::PeopleCreater(const string& file_name)
  : m_min_age(10),
    m_max_age(90),
    m_min_latitude(51.295449),
    m_max_latitude(51.709192),
    m_min_longtitude(-0.471357),
    m_max_longtitude(0.157610)
{
  m_first_names.reserve(50000);
  m_family_names.reserve(50000);

  ifstream in(file_name.c_str());
  if (!in.is_open()) 
    throw runtime_error("Cannot open " + file_name);

  typedef tokenizer< escaped_list_separator<char> > Tokenizer;

  vector< string > vec;
  string line;
  getline(in,line);
  while (getline(in,line))
    {
      Tokenizer tok(line);
      vec.assign(tok.begin(),tok.end());

      if (vec.size() != 2)
        throw runtime_error(line);
      
      m_first_names.push_back(string(vec[0]));
      m_family_names.push_back(string(vec[1]));
    }
}

boost::shared_ptr<KdPointConvertable> PeopleCreater::generate_a_person() const
{
  const string first_name = m_first_names[static_cast<double>(rand()) / static_cast<double>(RAND_MAX) * (m_first_names.size()-1)];
  const string family_name = m_family_names[static_cast<double>(rand()) / static_cast<double>(RAND_MAX) * (m_family_names.size()-1)];

  const static int age_range = m_max_age - m_min_age;
  const int age = static_cast<double>(rand()) / static_cast<double>(RAND_MAX) * age_range + m_min_age;

  const static double latitude_range = m_max_latitude - m_min_latitude;
  const double latitude = static_cast<double>(rand()) / static_cast<double>(RAND_MAX) * latitude_range + m_min_latitude;

  const static double longtitude_range = m_max_longtitude - m_min_longtitude;
  const double longtitude = static_cast<double>(rand()) / static_cast<double>(RAND_MAX) * longtitude_range + m_min_longtitude;
 
  return boost::shared_ptr<KdPointConvertable>(new Person(age, latitude, longtitude, string(first_name + " " + family_name))); 
}

vector<boost::shared_ptr<KdPointConvertable> > PeopleCreater::generate_people(const size_t number_of_people) const
{
  vector<boost::shared_ptr<KdPointConvertable> > people;
  people.reserve(number_of_people);
  for(int i = 0; i < number_of_people; i++)
    people.push_back(generate_a_person());
  return people;
}

void print_help()
{
  cout<<"Usage: NN_finer [data_size (in million)] [fast_launch (1 for yes)] [weight_for_age_difference] [weight_for_geography_distance]"<<endl;
  cout<<"Note: weight_for_age_difference and weight_for_geography_distance need to be either both provided or both omitted."<<endl;
}

bool is_input_valid(const int age, const double latitude, const double longtitude)
{
  if(age <= 0 || age > 130)
    return false;

  if(latitude > 180 || latitude <= -180)
    return false;

  if(longtitude > 180 || longtitude <= -180)
    return false;

  return true;
}

int main(int argc, char **argv)
{
  if(argc==2 && strcmp(argv[1], "help")==0){
    print_help();
    return 0;
  }
    
  const int data_size_million = argc > 1 ? atoi(argv[1]) : 1;
  const int fast_launch = argc > 2 ? atoi(argv[2]) : 1;
  if(fast_launch == 1){
    KdTree::set_get_partition_dimension_algorithm(KdTree::ROUND_ROBIN);
    KdTree::set_partition_algorithm(KdTree::RANDOM);
  }else{
    KdTree::set_get_partition_dimension_algorithm(KdTree::LARGEST_VARIANCE);
    KdTree::set_partition_algorithm(KdTree::QUICK_SELECT);
  }

  if(argc > 3){
    if(argc==3){
      print_help();
      return 0;
    }

    const double weight_age = argc > 3 ? atoi(argv[3]) : 1;
    const double weight_geography = argc > 4 ? atoi(argv[4]) : 10;
    vector<double> weights;
    weights.push_back(weight_age);
    weights.push_back(weight_geography);
    weights.push_back(weight_geography);
    Distance::set_weights(weights);
  }

  cout<<"Launching Nearest Neighbours Finder ..."<<endl<<flush;
  boost::shared_ptr<Normaliser> normaliser(new Normaliser());
  PeopleCreater people_generator;
  vector<shared_ptr<KdPointConvertable> > people(people_generator.generate_people(data_size_million * 1000000));
  KdTree kd_tree(people, normaliser);
  cout<<"Launching finished."<<endl;
  while(true){
    cout<<"Please enter your age, latitude, and longtitude"<<endl;
    int age = -1;
    double latitude = 181;
    double longtitude = 181;
    
    cin.clear();
    cin.ignore(numeric_limits<streamsize>::max(), '\n');
    
    cin>>age;
    if(cin.fail())
      continue;

    cin>>latitude;
    if(cin.fail())
      continue;

    cin>>longtitude;
    if(cin.fail())
      continue;
    
    if(!is_input_valid(age, latitude, longtitude)){
      cin.clear();
      cin.ignore(numeric_limits<streamsize>::max(), '\n');
      continue;
    }
    boost::shared_ptr<KdPointConvertable> person(new Person(age, latitude, longtitude, "Anonymous"));
    cout<<"Your info: ";
    person->print(cout);

    char c;
    const vector<boost::shared_ptr<KdPointConvertable> > results(kd_tree.get_NN(10, person.get()));
    cout<<"Your nearest neighbours are,"<<endl;
    for(int i = 0; i < results.size(); i++){
      results[i]->print(cout);
    }
    cout<<endl<<"Entering Q for quit, others for continue."<<endl;
    cin>>c;
    if(c == 'Q' || c == 'q')
      return 0;
  }                                      
}