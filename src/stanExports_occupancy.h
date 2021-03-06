// Generated by rstantools.  Do not edit by hand.

/*
    rstanOcc is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    rstanOcc is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with rstanOcc.  If not, see <http://www.gnu.org/licenses/>.
*/
#ifndef MODELS_HPP
#define MODELS_HPP
#define STAN__SERVICES__COMMAND_HPP
#include <rstan/rstaninc.hpp>
// Code generated by Stan version 2.19.1
#include <stan/model/model_header.hpp>
namespace model_occupancy_namespace {
using std::istream;
using std::string;
using std::stringstream;
using std::vector;
using stan::io::dump;
using stan::math::lgamma;
using stan::model::prob_grad;
using namespace stan::math;
static int current_statement_begin__;
stan::io::program_reader prog_reader__() {
    stan::io::program_reader reader;
    reader.add_event(0, 0, "start", "model_occupancy");
    reader.add_event(44, 42, "end", "model_occupancy");
    return reader;
}
#include <stan_meta_header.hpp>
class model_occupancy : public prob_grad {
private:
        int n_sampling_units;
        int n_surveys;
        std::vector<std::vector<int> > y;
public:
    model_occupancy(stan::io::var_context& context__,
        std::ostream* pstream__ = 0)
        : prob_grad(0) {
        ctor_body(context__, 0, pstream__);
    }
    model_occupancy(stan::io::var_context& context__,
        unsigned int random_seed__,
        std::ostream* pstream__ = 0)
        : prob_grad(0) {
        ctor_body(context__, random_seed__, pstream__);
    }
    void ctor_body(stan::io::var_context& context__,
                   unsigned int random_seed__,
                   std::ostream* pstream__) {
        typedef double local_scalar_t__;
        boost::ecuyer1988 base_rng__ =
          stan::services::util::create_rng(random_seed__, 0);
        (void) base_rng__;  // suppress unused var warning
        current_statement_begin__ = -1;
        static const char* function__ = "model_occupancy_namespace::model_occupancy";
        (void) function__;  // dummy to suppress unused var warning
        size_t pos__;
        (void) pos__;  // dummy to suppress unused var warning
        std::vector<int> vals_i__;
        std::vector<double> vals_r__;
        local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
        (void) DUMMY_VAR__;  // suppress unused var warning
        try {
            // initialize data block variables from context__
            current_statement_begin__ = 15;
            context__.validate_dims("data initialization", "n_sampling_units", "int", context__.to_vec());
            n_sampling_units = int(0);
            vals_i__ = context__.vals_i("n_sampling_units");
            pos__ = 0;
            n_sampling_units = vals_i__[pos__++];
            check_greater_or_equal(function__, "n_sampling_units", n_sampling_units, 0);
            current_statement_begin__ = 16;
            context__.validate_dims("data initialization", "n_surveys", "int", context__.to_vec());
            n_surveys = int(0);
            vals_i__ = context__.vals_i("n_surveys");
            pos__ = 0;
            n_surveys = vals_i__[pos__++];
            check_greater_or_equal(function__, "n_surveys", n_surveys, 0);
            current_statement_begin__ = 17;
            validate_non_negative_index("y", "n_sampling_units", n_sampling_units);
            validate_non_negative_index("y", "n_surveys", n_surveys);
            context__.validate_dims("data initialization", "y", "int", context__.to_vec(n_sampling_units,n_surveys));
            y = std::vector<std::vector<int> >(n_sampling_units, std::vector<int>(n_surveys, int(0)));
            vals_i__ = context__.vals_i("y");
            pos__ = 0;
            size_t y_k_0_max__ = n_sampling_units;
            size_t y_k_1_max__ = n_surveys;
            for (size_t k_1__ = 0; k_1__ < y_k_1_max__; ++k_1__) {
                for (size_t k_0__ = 0; k_0__ < y_k_0_max__; ++k_0__) {
                    y[k_0__][k_1__] = vals_i__[pos__++];
                }
            }
            size_t y_i_0_max__ = n_sampling_units;
            size_t y_i_1_max__ = n_surveys;
            for (size_t i_0__ = 0; i_0__ < y_i_0_max__; ++i_0__) {
                for (size_t i_1__ = 0; i_1__ < y_i_1_max__; ++i_1__) {
                    check_greater_or_equal(function__, "y[i_0__][i_1__]", y[i_0__][i_1__], 0);
                    check_less_or_equal(function__, "y[i_0__][i_1__]", y[i_0__][i_1__], 1);
                }
            }
            // initialize transformed data variables
            // execute transformed data statements
            // validate transformed data
            // validate, set parameter ranges
            num_params_r__ = 0U;
            param_ranges_i__.clear();
            current_statement_begin__ = 20;
            num_params_r__ += 1;
            current_statement_begin__ = 21;
            num_params_r__ += 1;
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(e, current_statement_begin__, prog_reader__());
            // Next line prevents compiler griping about no return
            throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***");
        }
    }
    ~model_occupancy() { }
    void transform_inits(const stan::io::var_context& context__,
                         std::vector<int>& params_i__,
                         std::vector<double>& params_r__,
                         std::ostream* pstream__) const {
        typedef double local_scalar_t__;
        stan::io::writer<double> writer__(params_r__, params_i__);
        size_t pos__;
        (void) pos__; // dummy call to supress warning
        std::vector<double> vals_r__;
        std::vector<int> vals_i__;
        current_statement_begin__ = 20;
        if (!(context__.contains_r("psi")))
            stan::lang::rethrow_located(std::runtime_error(std::string("Variable psi missing")), current_statement_begin__, prog_reader__());
        vals_r__ = context__.vals_r("psi");
        pos__ = 0U;
        context__.validate_dims("parameter initialization", "psi", "double", context__.to_vec());
        double psi(0);
        psi = vals_r__[pos__++];
        try {
            writer__.scalar_lub_unconstrain(0, 1, psi);
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(std::runtime_error(std::string("Error transforming variable psi: ") + e.what()), current_statement_begin__, prog_reader__());
        }
        current_statement_begin__ = 21;
        if (!(context__.contains_r("p")))
            stan::lang::rethrow_located(std::runtime_error(std::string("Variable p missing")), current_statement_begin__, prog_reader__());
        vals_r__ = context__.vals_r("p");
        pos__ = 0U;
        context__.validate_dims("parameter initialization", "p", "double", context__.to_vec());
        double p(0);
        p = vals_r__[pos__++];
        try {
            writer__.scalar_lub_unconstrain(0, 1, p);
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(std::runtime_error(std::string("Error transforming variable p: ") + e.what()), current_statement_begin__, prog_reader__());
        }
        params_r__ = writer__.data_r();
        params_i__ = writer__.data_i();
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
    template <bool propto__, bool jacobian__, typename T__>
    T__ log_prob(std::vector<T__>& params_r__,
                 std::vector<int>& params_i__,
                 std::ostream* pstream__ = 0) const {
        typedef T__ local_scalar_t__;
        local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
        (void) DUMMY_VAR__;  // dummy to suppress unused var warning
        T__ lp__(0.0);
        stan::math::accumulator<T__> lp_accum__;
        try {
            stan::io::reader<local_scalar_t__> in__(params_r__, params_i__);
            // model parameters
            current_statement_begin__ = 20;
            local_scalar_t__ psi;
            (void) psi;  // dummy to suppress unused var warning
            if (jacobian__)
                psi = in__.scalar_lub_constrain(0, 1, lp__);
            else
                psi = in__.scalar_lub_constrain(0, 1);
            current_statement_begin__ = 21;
            local_scalar_t__ p;
            (void) p;  // dummy to suppress unused var warning
            if (jacobian__)
                p = in__.scalar_lub_constrain(0, 1, lp__);
            else
                p = in__.scalar_lub_constrain(0, 1);
            // model body
            {
            current_statement_begin__ = 25;
            local_scalar_t__ log_psi(DUMMY_VAR__);
            (void) log_psi;  // dummy to suppress unused var warning
            stan::math::initialize(log_psi, DUMMY_VAR__);
            stan::math::fill(log_psi, DUMMY_VAR__);
            current_statement_begin__ = 26;
            local_scalar_t__ log1m_psi(DUMMY_VAR__);
            (void) log1m_psi;  // dummy to suppress unused var warning
            stan::math::initialize(log1m_psi, DUMMY_VAR__);
            stan::math::fill(log1m_psi, DUMMY_VAR__);
            current_statement_begin__ = 27;
            stan::math::assign(log_psi, stan::math::log(psi));
            current_statement_begin__ = 28;
            stan::math::assign(log1m_psi, log1m(psi));
            current_statement_begin__ = 31;
            lp_accum__.add(uniform_log<propto__>(psi, 0, 1));
            current_statement_begin__ = 32;
            lp_accum__.add(uniform_log<propto__>(p, 0, 1));
            current_statement_begin__ = 35;
            for (int r = 1; r <= n_sampling_units; ++r) {
                current_statement_begin__ = 36;
                if (as_bool(logical_gt(sum(get_base1(y, r, "y", 1)), 0))) {
                    current_statement_begin__ = 37;
                    lp_accum__.add((log_psi + bernoulli_log(get_base1(y, r, "y", 1), p)));
                } else {
                    current_statement_begin__ = 39;
                    lp_accum__.add(log_sum_exp((log_psi + bernoulli_log(get_base1(y, r, "y", 1), p)), log1m_psi));
                }
            }
            }
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(e, current_statement_begin__, prog_reader__());
            // Next line prevents compiler griping about no return
            throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***");
        }
        lp_accum__.add(lp__);
        return lp_accum__.sum();
    } // log_prob()
    template <bool propto, bool jacobian, typename T_>
    T_ log_prob(Eigen::Matrix<T_,Eigen::Dynamic,1>& params_r,
               std::ostream* pstream = 0) const {
      std::vector<T_> vec_params_r;
      vec_params_r.reserve(params_r.size());
      for (int i = 0; i < params_r.size(); ++i)
        vec_params_r.push_back(params_r(i));
      std::vector<int> vec_params_i;
      return log_prob<propto,jacobian,T_>(vec_params_r, vec_params_i, pstream);
    }
    void get_param_names(std::vector<std::string>& names__) const {
        names__.resize(0);
        names__.push_back("psi");
        names__.push_back("p");
    }
    void get_dims(std::vector<std::vector<size_t> >& dimss__) const {
        dimss__.resize(0);
        std::vector<size_t> dims__;
        dims__.resize(0);
        dimss__.push_back(dims__);
        dims__.resize(0);
        dimss__.push_back(dims__);
    }
    template <typename RNG>
    void write_array(RNG& base_rng__,
                     std::vector<double>& params_r__,
                     std::vector<int>& params_i__,
                     std::vector<double>& vars__,
                     bool include_tparams__ = true,
                     bool include_gqs__ = true,
                     std::ostream* pstream__ = 0) const {
        typedef double local_scalar_t__;
        vars__.resize(0);
        stan::io::reader<local_scalar_t__> in__(params_r__, params_i__);
        static const char* function__ = "model_occupancy_namespace::write_array";
        (void) function__;  // dummy to suppress unused var warning
        // read-transform, write parameters
        double psi = in__.scalar_lub_constrain(0, 1);
        vars__.push_back(psi);
        double p = in__.scalar_lub_constrain(0, 1);
        vars__.push_back(p);
        double lp__ = 0.0;
        (void) lp__;  // dummy to suppress unused var warning
        stan::math::accumulator<double> lp_accum__;
        local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
        (void) DUMMY_VAR__;  // suppress unused var warning
        if (!include_tparams__ && !include_gqs__) return;
        try {
            if (!include_gqs__ && !include_tparams__) return;
            if (!include_gqs__) return;
        } catch (const std::exception& e) {
            stan::lang::rethrow_located(e, current_statement_begin__, prog_reader__());
            // Next line prevents compiler griping about no return
            throw std::runtime_error("*** IF YOU SEE THIS, PLEASE REPORT A BUG ***");
        }
    }
    template <typename RNG>
    void write_array(RNG& base_rng,
                     Eigen::Matrix<double,Eigen::Dynamic,1>& params_r,
                     Eigen::Matrix<double,Eigen::Dynamic,1>& vars,
                     bool include_tparams = true,
                     bool include_gqs = true,
                     std::ostream* pstream = 0) const {
      std::vector<double> params_r_vec(params_r.size());
      for (int i = 0; i < params_r.size(); ++i)
        params_r_vec[i] = params_r(i);
      std::vector<double> vars_vec;
      std::vector<int> params_i_vec;
      write_array(base_rng, params_r_vec, params_i_vec, vars_vec, include_tparams, include_gqs, pstream);
      vars.resize(vars_vec.size());
      for (int i = 0; i < vars.size(); ++i)
        vars(i) = vars_vec[i];
    }
    static std::string model_name() {
        return "model_occupancy";
    }
    void constrained_param_names(std::vector<std::string>& param_names__,
                                 bool include_tparams__ = true,
                                 bool include_gqs__ = true) const {
        std::stringstream param_name_stream__;
        param_name_stream__.str(std::string());
        param_name_stream__ << "psi";
        param_names__.push_back(param_name_stream__.str());
        param_name_stream__.str(std::string());
        param_name_stream__ << "p";
        param_names__.push_back(param_name_stream__.str());
        if (!include_gqs__ && !include_tparams__) return;
        if (include_tparams__) {
        }
        if (!include_gqs__) return;
    }
    void unconstrained_param_names(std::vector<std::string>& param_names__,
                                   bool include_tparams__ = true,
                                   bool include_gqs__ = true) const {
        std::stringstream param_name_stream__;
        param_name_stream__.str(std::string());
        param_name_stream__ << "psi";
        param_names__.push_back(param_name_stream__.str());
        param_name_stream__.str(std::string());
        param_name_stream__ << "p";
        param_names__.push_back(param_name_stream__.str());
        if (!include_gqs__ && !include_tparams__) return;
        if (include_tparams__) {
        }
        if (!include_gqs__) return;
    }
}; // model
}  // namespace
typedef model_occupancy_namespace::model_occupancy stan_model;
#endif
