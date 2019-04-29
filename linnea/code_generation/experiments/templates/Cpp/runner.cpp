#include <array>
#include <string>
#include <utility>
#include <fstream>

#include <benchmarker/benchmarker.hpp>
#include <libraries.hpp>
#include <generator/generator.hpp>

#include "reference/naive_armadillo.hpp"
#include "reference/recommended_armadillo.hpp"
#include "reference/naive_eigen.hpp"
#include "reference/recommended_eigen.hpp"

#include "operand_generator.hpp"

typedef double float_type;

template<typename Function, typename Tuple, std::size_t... I>
decltype(auto) call(Function && f, Tuple && t, std::index_sequence<I...>)
{{
    return f(std::get<I>(t)...);
}}

int main(int argc, char ** argv)
{{
    std::cout << "Test runner for algorithm {0}" << std::endl;
    linalg_tests::basic_benchmarker<std::chrono::duration<float>> benchmark;
    benchmark.set_cache_size(30 * 1024 * 1024);
    std::ofstream file("cpp_results_{0}_timings.txt");
    {{
    generator::generator<library::arma, float_type> arma_gen{{100}};
    auto matrices = operand_generator(arma_gen);
    constexpr std::size_t tuple_length = std::tuple_size<decltype(matrices)>::value;
    typedef std::make_index_sequence<tuple_length> Indices;
    auto res_naive = call(naive_armadillo{{}}, matrices, Indices{{}});
    auto res_recomm = call(recommended_armadillo{{}}, matrices, Indices{{}});
    //std::cout << "Armadillo norm(naive-recom): " << arma::norm(res_naive - res_recomm) << std::endl;
    //std::cout << "Armadillo naive(0,0): " <<res_naive(0, 0) << std::endl;

    std::array<std::string, 1> labels_arma_n{{ {{"naive_armadillo"}} }};
    benchmark.run(file, labels_arma_n, 20, [&]() {{
            call(naive_armadillo{{}}, matrices, Indices{{}});
            }});

    std::array<std::string, 1> labels_arma_r{{ {{"recommended_armadillo"}} }};
    benchmark.run(file, labels_arma_r, 20, [&]() {{
            call(recommended_armadillo{{}}, matrices, Indices{{}});
            }});
    }}
    {{ 
    generator::generator<library::eigen, float_type> eigen_gen{{100}};
    auto matrices = operand_generator(eigen_gen);
    constexpr std::size_t tuple_length = std::tuple_size<decltype(matrices)>::value;
    typedef std::make_index_sequence<tuple_length> Indices;
    auto res_naive = call(naive_eigen{{}}, matrices, Indices{{}});
    auto res_recomm = call(recommended_eigen{{}}, matrices, Indices{{}});
    //std::cout << "Eigen norm(naive-recom): " << (res_naive - res_recomm).norm() << std::endl;
    //std::cout << "Eigen naive(0,0): " << res_naive(0, 0) << std::endl;

    std::array<std::string, 1> labels_eigen_n{{ {{"naive_eigen"}} }};
    benchmark.run(file, labels_eigen_n, 20, [&]() {{
            call(naive_eigen{{}}, matrices, Indices{{}});
            }});

    std::array<std::string, 1> labels_eigen_r{{ {{"recommended_eigen"}} }};
    benchmark.run(file, labels_eigen_r, 20, [&]() {{
            call(recommended_eigen{{}}, matrices, Indices{{}});
            }});
    }}
    
    file.close();
}}