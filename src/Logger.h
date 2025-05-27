#ifndef LOGGER_H
#define LOGGER_H

#include <iostream>
#include <fstream>
#include <string>
#include <chrono>
#include <iomanip>
#include <mutex>
#include <memory>
#include <sstream>
#include <format>

namespace libqm {
class Logger {
public:
  enum class LogLevel {
    DEBUG,
    INFO,
    WARNING,
    ERROR
  };

  Logger(const Logger&) = delete;
  Logger& operator=(const Logger&) = delete;
  Logger(Logger&&) = delete;
  Logger& operator=(Logger&&) = delete;

  static Logger& the() {
    static Logger instance;
    return instance;
  }

  void set_level(LogLevel level);
  void set_to_console(bool enable);
  bool set_file(const std::string& filename);
  void log(LogLevel level, const std::string& message);
  void debug(const std::string& message);
  void info(const std::string& message);
  void warn(const std::string& message);
  void error(const std::string& message);

  template<typename... Args>
  void log_fmt(LogLevel level, std::format_string<Args...> fmt, Args&&... args) {
    try {
      std::string formatted_message = std::format(fmt, std::forward<Args>(args)...);
      log(level, formatted_message);
    } catch (const std::format_error& e) {
      log(LogLevel::ERROR, std::string("Formatting error: ") + e.what());
    }
  }

  template<typename... Args> void debug_fmt(std::format_string<Args...> fmt, Args&&... args) { log_fmt(LogLevel::DEBUG, fmt, std::forward<Args>(args)...); }
  template<typename... Args> void info_fmt(std::format_string<Args...> fmt, Args&&... args)  { log_fmt(LogLevel::INFO, fmt, std::forward<Args>(args)...); }
  template<typename... Args> void warn_fmt(std::format_string<Args...> fmt, Args&&... args)  { log_fmt(LogLevel::WARNING, fmt, std::forward<Args>(args)...); }
  template<typename... Args> void error_fmt(std::format_string<Args...> fmt, Args&&... args) { logFmt(LogLevel::ERROR, fmt, std::forward<Args>(args)...); }

  template<typename T, typename... Args>
  void log_stream(LogLevel level, const T& firstArg, const Args&... restArgs) {
    std::ostringstream oss;
    oss << firstArg;
    ((oss << " " << restArgs), ...);
    log(level, oss.str());
  }

  template<typename T, typename... Args> void debug_stream(const T& firstArg, const Args&... restArgs) { log_stream(LogLevel::DEBUG, firstArg, restArgs...); }
  template<typename T, typename... Args> void info_stream(const T& firstArg, const Args&... restArgs)  { log_stream(LogLevel::INFO, firstArg, restArgs...); }
  template<typename T, typename... Args> void warn_stream(const T& firstArg, const Args&... restArgs)  { log_stream(LogLevel::WARNING, firstArg, restArgs...); }
  template<typename T, typename... Args> void error_stream(const T& firstArg, const Args&... restArgs) { logStream(LogLevel::ERROR, firstArg, restArgs...); }

private:
  Logger();
  ~Logger() = default;

  std::string level_to_string(LogLevel level) const;

  LogLevel level_;
  std::unique_ptr<std::ofstream> file_;
  std::string path_;
  bool to_console_;
  std::mutex mutex_;
};
} // namespace libqm

#endif // LOGGER_H
