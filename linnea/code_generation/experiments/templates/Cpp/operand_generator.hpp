#include <generator/generator.hpp>

template<typename Gen>
decltype(auto) operand_generator(Gen && gen)
{{
{}
    return std::make_tuple({});
}}