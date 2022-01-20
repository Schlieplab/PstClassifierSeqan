#pragma once

#include <cmath>
#include <functional>
#include <iostream>
#include <thread>
#include <vector>

namespace pst::parallelize {

enum ExecutionStrategy { parallel, sequential };

std::vector<std::tuple<size_t, size_t>> get_bounds(size_t size) {
  const size_t processor_count = std::thread::hardware_concurrency();
  std::vector<std::tuple<size_t, size_t>> bounds_per_thread{};
  float values_per_thread = float(size) / float(processor_count);

  for (size_t i = 0; i < std::min(processor_count, size); i++) {
    size_t start_index = std::floor(values_per_thread * i);
    size_t stop_index = std::ceil(values_per_thread * (i + 1.0));
    if (i == (processor_count - 1)) {
      stop_index = size;
    }
    bounds_per_thread.emplace_back(start_index, stop_index);
  }

  return bounds_per_thread;
}

void parallelize(size_t size, const std::function<void(size_t, size_t)> &fun) {
  std::vector<std::thread> threads{};

  auto bounds = get_bounds(size);
  for (auto &[start_index, stop_index] : bounds) {
    threads.emplace_back(fun, start_index, stop_index);
  }

  for (auto &thread : threads) {
    if (thread.joinable()) {
      thread.join();
    }
  }
}

void parallelize_with_progress(
    size_t size,
    const std::function<void(
        size_t, size_t, indicators::DynamicProgress<indicators::ProgressBar> &,
        int)> &fun) {
  std::vector<std::thread> threads{};

  indicators::DynamicProgress<indicators::ProgressBar> bars{};
  bars.set_option(indicators::option::HideBarWhenComplete{true});

  indicators::show_console_cursor(false);

  auto bounds = get_bounds(size);
  for (auto &[start_index, stop_index] : bounds) {
    indicators::ProgressBar bar{indicators::option::BarWidth{50},
                                indicators::option::PrefixText{"Distances: "},
                                indicators::option::ShowElapsedTime{true},
                                indicators::option::ShowRemainingTime{true}};
    auto bar_i = bars.push_back(bar);
    threads.emplace_back(fun, start_index, stop_index, std::ref(bars), bar_i);
  }

  for (auto &thread : threads) {
    if (thread.joinable()) {
      thread.join();
    }
  }
  indicators::show_console_cursor(true);
}

} // namespace pst::parallelize
