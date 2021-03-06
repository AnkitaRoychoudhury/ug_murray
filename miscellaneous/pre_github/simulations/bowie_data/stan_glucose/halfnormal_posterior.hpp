
// Code generated by stanc f556d0d
#include <stan/model/model_header.hpp>
namespace halfnormal_posterior_model_namespace {

template <typename T, typename S>
std::vector<T> resize_to_match__(std::vector<T>& dst, const std::vector<S>& src) {
  dst.resize(src.size());
  return dst;
}

template <typename T>
Eigen::Matrix<T, -1, -1>
resize_to_match__(Eigen::Matrix<T, -1, -1>& dst, const Eigen::Matrix<T, -1, -1>& src) {
  dst.resize(src.rows(), src.cols());
  return dst;
}

template <typename T>
Eigen::Matrix<T, 1, -1>
resize_to_match__(Eigen::Matrix<T, 1, -1>& dst, const Eigen::Matrix<T, 1, -1>& src) {
  dst.resize(src.size());
  return dst;
}

template <typename T>
Eigen::Matrix<T, -1, 1>
resize_to_match__(Eigen::Matrix<T, -1, 1>& dst, const Eigen::Matrix<T, -1, 1>& src) {
  dst.resize(src.size());
  return dst;
}
std::vector<double> to_doubles__(std::initializer_list<double> x) {
  return x;
}

std::vector<stan::math::var> to_vars__(std::initializer_list<stan::math::var> x) {
  return x;
}

inline void validate_positive_index(const char* var_name, const char* expr,
                                    int val) {
  if (val < 1) {
    std::stringstream msg;
    msg << "Found dimension size less than one in simplex declaration"
        << "; variable=" << var_name << "; dimension size expression=" << expr
        << "; expression value=" << val;
    std::string msg_str(msg.str());
    throw std::invalid_argument(msg_str.c_str());
  }
}

inline void validate_unit_vector_index(const char* var_name, const char* expr,
                                       int val) {
  if (val <= 1) {
    std::stringstream msg;
    if (val == 1) {
      msg << "Found dimension size one in unit vector declaration."
          << " One-dimensional unit vector is discrete"
          << " but the target distribution must be continuous."
          << " variable=" << var_name << "; dimension size expression=" << expr;
    } else {
      msg << "Found dimension size less than one in unit vector declaration"
          << "; variable=" << var_name << "; dimension size expression=" << expr
          << "; expression value=" << val;
    }
    std::string msg_str(msg.str());
    throw std::invalid_argument(msg_str.c_str());
  }
}


using std::istream;
using std::string;
using std::stringstream;
using std::vector;
using std::pow;
using stan::io::dump;
using stan::math::lgamma;
using stan::model::model_base_crtp;
using stan::model::rvalue;
using stan::model::cons_list;
using stan::model::index_uni;
using stan::model::index_max;
using stan::model::index_min;
using stan::model::index_min_max;
using stan::model::index_multi;
using stan::model::index_omni;
using stan::model::nil_index_list;
using namespace stan::math; 

static int current_statement__ = 0;
static const std::vector<string> locations_array__ = {" (found before start of program)",
                                                      " (in '/Users/ankitaroychoudhury/Documents/MURRAY/simulations/bowie_data/stan_glucose/halfnormal_posterior.stan', line 8, column 2 to column 19)",
                                                      " (in '/Users/ankitaroychoudhury/Documents/MURRAY/simulations/bowie_data/stan_glucose/halfnormal_posterior.stan', line 9, column 2 to column 22)",
                                                      " (in '/Users/ankitaroychoudhury/Documents/MURRAY/simulations/bowie_data/stan_glucose/halfnormal_posterior.stan', line 27, column 4 to column 22)",
                                                      " (in '/Users/ankitaroychoudhury/Documents/MURRAY/simulations/bowie_data/stan_glucose/halfnormal_posterior.stan', line 30, column 8 to column 50)",
                                                      " (in '/Users/ankitaroychoudhury/Documents/MURRAY/simulations/bowie_data/stan_glucose/halfnormal_posterior.stan', line 29, column 22 to line 33, column 9)",
                                                      " (in '/Users/ankitaroychoudhury/Documents/MURRAY/simulations/bowie_data/stan_glucose/halfnormal_posterior.stan', line 29, column 4 to line 33, column 9)",
                                                      " (in '/Users/ankitaroychoudhury/Documents/MURRAY/simulations/bowie_data/stan_glucose/halfnormal_posterior.stan', line 15, column 2 to column 24)",
                                                      " (in '/Users/ankitaroychoudhury/Documents/MURRAY/simulations/bowie_data/stan_glucose/halfnormal_posterior.stan', line 16, column 2 to column 26)",
                                                      " (in '/Users/ankitaroychoudhury/Documents/MURRAY/simulations/bowie_data/stan_glucose/halfnormal_posterior.stan', line 22, column 2 to column 24)",
                                                      " (in '/Users/ankitaroychoudhury/Documents/MURRAY/simulations/bowie_data/stan_glucose/halfnormal_posterior.stan', line 2, column 2 to column 17)",
                                                      " (in '/Users/ankitaroychoudhury/Documents/MURRAY/simulations/bowie_data/stan_glucose/halfnormal_posterior.stan', line 3, column 2 to column 12)",
                                                      " (in '/Users/ankitaroychoudhury/Documents/MURRAY/simulations/bowie_data/stan_glucose/halfnormal_posterior.stan', line 4, column 2 to column 21)"};



class halfnormal_posterior_model : public model_base_crtp<halfnormal_posterior_model> {

 private:
  int pos__;
  int N;
  std::vector<double> k;
  int N_ppc;
 
 public:
  ~halfnormal_posterior_model() { }
  
  std::string model_name() const { return "halfnormal_posterior_model"; }
  
  halfnormal_posterior_model(stan::io::var_context& context__,
                             unsigned int random_seed__ = 0,
                             std::ostream* pstream__ = nullptr) : model_base_crtp(0) {
    typedef double local_scalar_t__;
    boost::ecuyer1988 base_rng__ = 
        stan::services::util::create_rng(random_seed__, 0);
    (void) base_rng__;  // suppress unused var warning
    static const char* function__ = "halfnormal_posterior_model_namespace::halfnormal_posterior_model";
    (void) function__;  // suppress unused var warning
    
    try {
      
      pos__ = 1;
      context__.validate_dims("data initialization","N","int",
          context__.to_vec());
      
      current_statement__ = 10;
      N = context__.vals_i("N")[(1 - 1)];
      current_statement__ = 11;
      validate_non_negative_index("k", "N", N);
      context__.validate_dims("data initialization","k","double",
          context__.to_vec(N));
      k = std::vector<double>(N, 0);
      
      current_statement__ = 11;
      assign(k, nil_index_list(), context__.vals_r("k"),
        "assigning variable k");
      context__.validate_dims("data initialization","N_ppc","int",
          context__.to_vec());
      
      current_statement__ = 12;
      N_ppc = context__.vals_i("N_ppc")[(1 - 1)];
      current_statement__ = 10;
      current_statement__ = 10;
      check_greater_or_equal(function__, "N", N, 0);
      current_statement__ = 12;
      current_statement__ = 12;
      check_greater_or_equal(function__, "N_ppc", N_ppc, 0);
    } catch (const std::exception& e) {
      stan::lang::rethrow_located(e, locations_array__[current_statement__]);
      // Next line prevents compiler griping about no return
      throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***"); 
    }
    num_params_r__ = 0U;
    
    try {
      num_params_r__ += 1;
      num_params_r__ += 1;
    } catch (const std::exception& e) {
      stan::lang::rethrow_located(e, locations_array__[current_statement__]);
      // Next line prevents compiler griping about no return
      throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***"); 
    }
  }
  template <bool propto__, bool jacobian__, typename T__>
  T__ log_prob(std::vector<T__>& params_r__, std::vector<int>& params_i__,
               std::ostream* pstream__ = 0) const {
    typedef T__ local_scalar_t__;
    T__ lp__(0.0);
    stan::math::accumulator<T__> lp_accum__;
    static const char* function__ = "halfnormal_posterior_model_namespace::log_prob";
(void) function__;  // suppress unused var warning

    stan::io::reader<local_scalar_t__> in__(params_r__, params_i__);
    
    try {
      local_scalar_t__ mu;
      
      current_statement__ = 1;
      mu = in__.scalar();
      current_statement__ = 1;
      if (jacobian__) {
        current_statement__ = 1;
        mu = stan::math::lb_constrain(mu, 0, lp__);
      } else {
        current_statement__ = 1;
        mu = stan::math::lb_constrain(mu, 0);
      }
      local_scalar_t__ sigma;
      
      current_statement__ = 2;
      sigma = in__.scalar();
      current_statement__ = 2;
      if (jacobian__) {
        current_statement__ = 2;
        sigma = stan::math::lb_constrain(sigma, 0, lp__);
      } else {
        current_statement__ = 2;
        sigma = stan::math::lb_constrain(sigma, 0);
      }
      {
        current_statement__ = 7;
        lp_accum__.add(lognormal_log<propto__>(mu, 0, 20));
        current_statement__ = 8;
        lp_accum__.add(lognormal_log<propto__>(sigma, 0, 30));
        current_statement__ = 9;
        lp_accum__.add(normal_log<propto__>(k, mu, sigma));
      }
    } catch (const std::exception& e) {
      stan::lang::rethrow_located(e, locations_array__[current_statement__]);
      // Next line prevents compiler griping about no return
      throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***"); 
    }
    lp_accum__.add(lp__);
    return lp_accum__.sum();
    } // log_prob() 
    
  template <typename RNG>
  void write_array(RNG& base_rng__, std::vector<double>& params_r__,
                   std::vector<int>& params_i__, std::vector<double>& vars__,
                   bool emit_transformed_parameters__ = true,
                   bool emit_generated_quantities__ = true,
                   std::ostream* pstream__ = 0) const {
    typedef double local_scalar_t__;
    vars__.resize(0);
    stan::io::reader<local_scalar_t__> in__(params_r__, params_i__);
    static const char* function__ = "halfnormal_posterior_model_namespace::write_array";
(void) function__;  // suppress unused var warning

    (void) function__;  // suppress unused var warning

    double lp__ = 0.0;
    (void) lp__;  // dummy to suppress unused var warning
    stan::math::accumulator<double> lp_accum__;
    
    try {
      double mu;
      
      current_statement__ = 1;
      mu = in__.scalar();
      current_statement__ = 1;
      mu = stan::math::lb_constrain(mu, 0);
      double sigma;
      
      current_statement__ = 2;
      sigma = in__.scalar();
      current_statement__ = 2;
      sigma = stan::math::lb_constrain(sigma, 0);
      vars__.push_back(mu);
      vars__.push_back(sigma);
      if (logical_negation((primitive_value(emit_transformed_parameters__) ||
            primitive_value(emit_generated_quantities__)))) {
        return ;
      } 
      if (logical_negation(emit_generated_quantities__)) {
        return ;
      } 
      current_statement__ = 3;
      validate_non_negative_index("k_ppc", "N_ppc", N_ppc);
      std::vector<double> k_ppc;
      k_ppc = std::vector<double>(N_ppc, 0);
      
      current_statement__ = 6;
      for (size_t i = 1; i <= N_ppc; ++i) {
        current_statement__ = 4;
        assign(k_ppc, cons_list(index_uni(i), nil_index_list()),
          (mu + stan::math::abs(normal_rng(0, sigma, base_rng__))),
          "assigning variable k_ppc");}
      for (size_t sym1__ = 1; sym1__ <= N_ppc; ++sym1__) {
        vars__.push_back(k_ppc[(sym1__ - 1)]);}
    } catch (const std::exception& e) {
      stan::lang::rethrow_located(e, locations_array__[current_statement__]);
      // Next line prevents compiler griping about no return
      throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***"); 
    }
    } // write_array() 
    
  void transform_inits(const stan::io::var_context& context__,
                       std::vector<int>& params_i__,
                       std::vector<double>& vars__, std::ostream* pstream__) const {
    typedef double local_scalar_t__;
    vars__.resize(0);
    vars__.reserve(num_params_r__);
    
    try {
      int pos__;
      
      pos__ = 1;
      double mu;
      
      current_statement__ = 1;
      mu = context__.vals_r("mu")[(1 - 1)];
      current_statement__ = 1;
      mu = stan::math::lb_free(mu, 0);
      double sigma;
      
      current_statement__ = 2;
      sigma = context__.vals_r("sigma")[(1 - 1)];
      current_statement__ = 2;
      sigma = stan::math::lb_free(sigma, 0);
      vars__.push_back(mu);
      vars__.push_back(sigma);
    } catch (const std::exception& e) {
      stan::lang::rethrow_located(e, locations_array__[current_statement__]);
      // Next line prevents compiler griping about no return
      throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***"); 
    }
    } // transform_inits() 
    
  void get_param_names(std::vector<std::string>& names__) const {
    
    names__.resize(0);
    names__.push_back("mu");
    names__.push_back("sigma");
    names__.push_back("k_ppc");
    } // get_param_names() 
    
  void get_dims(std::vector<std::vector<size_t>>& dimss__) const {
    dimss__.resize(0);
    std::vector<size_t> dims__;
    dimss__.push_back(dims__);
    dims__.resize(0);
    dimss__.push_back(dims__);
    dims__.resize(0);
    dims__.push_back(N_ppc);
    dimss__.push_back(dims__);
    dims__.resize(0);
    
    } // get_dims() 
    
  void constrained_param_names(std::vector<std::string>& param_names__,
                               bool emit_transformed_parameters__ = true,
                               bool emit_generated_quantities__ = true) const {
    
    param_names__.push_back(std::string() + "mu");
    param_names__.push_back(std::string() + "sigma");
    if (emit_transformed_parameters__) {
      
    }
    
    if (emit_generated_quantities__) {
      for (size_t sym1__ = 1; sym1__ <= N_ppc; ++sym1__) {
        {
          param_names__.push_back(std::string() + "k_ppc" + '.' + std::to_string(sym1__));
        }}
    }
    
    } // constrained_param_names() 
    
  void unconstrained_param_names(std::vector<std::string>& param_names__,
                                 bool emit_transformed_parameters__ = true,
                                 bool emit_generated_quantities__ = true) const {
    
    param_names__.push_back(std::string() + "mu");
    param_names__.push_back(std::string() + "sigma");
    if (emit_transformed_parameters__) {
      
    }
    
    if (emit_generated_quantities__) {
      for (size_t sym1__ = 1; sym1__ <= N_ppc; ++sym1__) {
        {
          param_names__.push_back(std::string() + "k_ppc" + '.' + std::to_string(sym1__));
        }}
    }
    
    } // unconstrained_param_names() 
    
  std::string get_constrained_sizedtypes() const {
    stringstream s__;
    s__ << "[{\"name\":\"mu\",\"type\":{\"name\":\"real\"},\"block\":\"parameters\"},{\"name\":\"sigma\",\"type\":{\"name\":\"real\"},\"block\":\"parameters\"},{\"name\":\"k_ppc\",\"type\":{\"name\":\"array\",\"length\":" << N_ppc << ",\"element_type\":{\"name\":\"real\"}},\"block\":\"generated_quantities\"}]";
    return s__.str();
    } // get_constrained_sizedtypes() 
    
  std::string get_unconstrained_sizedtypes() const {
    stringstream s__;
    s__ << "[{\"name\":\"mu\",\"type\":{\"name\":\"real\"},\"block\":\"parameters\"},{\"name\":\"sigma\",\"type\":{\"name\":\"real\"},\"block\":\"parameters\"},{\"name\":\"k_ppc\",\"type\":{\"name\":\"array\",\"length\":" << N_ppc << ",\"element_type\":{\"name\":\"real\"}},\"block\":\"generated_quantities\"}]";
    return s__.str();
    } // get_unconstrained_sizedtypes() 
    
  
    // Begin method overload boilerplate
    template <typename RNG>
    void write_array(RNG& base_rng__,
                     Eigen::Matrix<double,Eigen::Dynamic,1>& params_r,
                     Eigen::Matrix<double,Eigen::Dynamic,1>& vars,
                     bool emit_transformed_parameters__ = true,
                     bool emit_generated_quantities__ = true,
                     std::ostream* pstream = 0) const {
      std::vector<double> params_r_vec(params_r.size());
      for (int i = 0; i < params_r.size(); ++i)
        params_r_vec[i] = params_r(i);
      std::vector<double> vars_vec;
      std::vector<int> params_i_vec;
      write_array(base_rng__, params_r_vec, params_i_vec, vars_vec,
          emit_transformed_parameters__, emit_generated_quantities__, pstream);
      vars.resize(vars_vec.size());
      for (int i = 0; i < vars.size(); ++i)
        vars(i) = vars_vec[i];
    }

    template <bool propto__, bool jacobian__, typename T_>
    T_ log_prob(Eigen::Matrix<T_,Eigen::Dynamic,1>& params_r,
               std::ostream* pstream = 0) const {
      std::vector<T_> vec_params_r;
      vec_params_r.reserve(params_r.size());
      for (int i = 0; i < params_r.size(); ++i)
        vec_params_r.push_back(params_r(i));
      std::vector<int> vec_params_i;
      return log_prob<propto__,jacobian__,T_>(vec_params_r, vec_params_i, pstream);
    }

    void transform_inits(const stan::io::var_context& context,
                         Eigen::Matrix<double, Eigen::Dynamic, 1>& params_r,
                         std::ostream* pstream__) const {
      std::vector<double> params_r_vec;
      std::vector<int> params_i_vec;
      transform_inits(context, params_i_vec, params_r_vec, pstream__);
      params_r.resize(params_r_vec.size());
      for (int i = 0; i < params_r.size(); ++i)
        params_r(i) = params_r_vec[i];
    }

};
}

typedef halfnormal_posterior_model_namespace::halfnormal_posterior_model stan_model;

#ifndef USING_R

// Boilerplate
stan::model::model_base& new_model(
        stan::io::var_context& data_context,
        unsigned int seed,
        std::ostream* msg_stream) {
  stan_model* m = new stan_model(data_context, seed, msg_stream);
  return *m;
}

#endif


