#pragma once

#include <list>
#include <utility>
#include <cstddef>
#include <optional>

template <
  class Key,
  class Value,
  template <typename...> class MapTypeTemplate
  >
class LruCache {
public:
  using key_type = Key;
  using value_type = Value;
  using list_type = std::list<key_type>;
  using list_iterator_type = typename list_type::iterator;
  using internal_value_type = std::pair<value_type, list_iterator_type>;
  using map_type = MapTypeTemplate<key_type, internal_value_type>;
  using map_iterator_type = typename map_type::iterator;
  using const_map_iterator_type = typename map_type::const_iterator;


  explicit LruCache(size_t capacity) : m_capacity(capacity) {}

  size_t size() const {
    return m_map.size();
  }

  size_t capacity() const {
    return m_capacity;
  }

  bool empty() const {
    return m_map.empty();
  }

  bool contains(const key_type &key) const {
    return m_map.find(key) != m_map.end();
  }

  void insert(const key_type &key, const value_type &value) {
    if (m_capacity == 0) return;

    map_iterator_type map_iter = m_map.find(key);

    if (map_iter == m_map.end()) {
      if (size() >= m_capacity) {
        evict();
      }
      m_list.push_front(key);
      m_map.emplace(key, std::make_pair(value, m_list.begin()));
    } else {
      map_iter->second.first = value;
      move_to_front(map_iter->second.second);
    }
  }

  std::optional<value_type> get(const key_type &key) {
    map_iterator_type map_iter = m_map.find(key);
    if (map_iter == m_map.end()) {
      return std::nullopt;
    }
    move_to_front(map_iter->second.second);
    return map_iter->second.first;
  }

  void clear() {
    m_map.clear();
    m_list.clear();
  }

private:
  void move_to_front(list_iterator_type list_iter) {
    if (list_iter != m_list.begin()) {
      m_list.splice(m_list.begin(), m_list, list_iter);
      map_iterator_type map_iter = m_map.find(*list_iter);
      if(map_iter != m_map.end()) {
        map_iter->second.second = m_list.begin();
      }
    }
  }


  void evict() {
    if (m_list.empty()) {
      return;
    }

    key_type key_to_evict = m_list.back();
    m_map.erase(key_to_evict);
    m_list.pop_back();
  }

private:
  map_type m_map;
  list_type m_list;
  size_t m_capacity;
};
