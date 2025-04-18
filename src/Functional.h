#pragma once

#include <type_traits>

#include "Expression.h"
#include "Operator.h"
#include "Term.h"

namespace libqm {
template <typename F, typename... Args>
using is_operator_callable = std::is_invocable<F, Operator, Args...>;

template <typename F, typename... Args>
auto transform_operator(F&& f, Operator op, Args&&... args)
    -> std::enable_if_t<is_operator_callable<F, Args...>::value,
                        decltype(f(op, std::forward<Args>(args)...))> {
  return f(op, std::forward<Args>(args)...);
}

template <typename F, typename... Args>
auto transform_term(F&& f, Term::complex_type c, const Term::container_type& ops, Args&&... args)
    -> std::enable_if_t<is_operator_callable<F, Args...>::value, Expression> {
  Expression result(c);
  for (const auto& op : ops) {
    result *= transform_operator(std::forward<F>(f), op, std::forward<Args>(args)...);
  }
  return result;
}

template <typename F, typename... Args>
auto transform_term(F&& f, Term::complex_type c, const Term& term, Args&&... args)
    -> std::enable_if_t<is_operator_callable<F, Args...>::value, Expression> {
  return transform_term(std::forward<F>(f), c, term.operators, std::forward<Args>(args)...);
}

template <typename F, typename... Args>
auto transform_expression(F&& f, const Expression& expr, Args&&... args)
    -> std::enable_if_t<is_operator_callable<F, Args...>::value, Expression> {
  Expression result;
  for (const auto& [ops, coeff] : expr.hashmap) {
    result += transform_term(std::forward<F>(f), coeff, ops, std::forward<Args>(args)...);
  }
  return result;
}
}  // namespace libqm