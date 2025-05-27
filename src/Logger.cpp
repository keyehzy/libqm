#include "Logger.h"
#include <stdexcept>
#include <ctime>

namespace libqm {
Logger::Logger() : level_(LogLevel::INFO), to_console_(true), file_(nullptr) {
    set_file("application.log");
}

void Logger::set_level(LogLevel level) {
    std::lock_guard<std::mutex> lock(mutex_);
    level_ = level;
}

void Logger::set_to_console(bool enable) {
    std::lock_guard<std::mutex> lock(mutex_);
    to_console_ = enable;
}

bool Logger::set_file(const std::string& filename) {
  std::lock_guard<std::mutex> lock(mutex_);
  path_ = filename;
  auto newFile = std::make_unique<std::ofstream>(path_, std::ios::app);

  if (newFile && newFile->is_open()) {
    file_ = std::move(newFile);
    return true;
  } else {
    std::cerr << "Error: Could not open log file: " << path_ << std::endl;
    file_.reset();
    return false;
  }
}

std::string Logger::level_to_string(LogLevel level) const {
  switch (level) {
  case LogLevel::DEBUG:   return "DEBUG";
  case LogLevel::INFO:    return "INFO";
  case LogLevel::WARNING: return "WARN";
  case LogLevel::ERROR:   return "ERROR";
  default:                return "UNKNOWN";
  }
}

void Logger::log(LogLevel level, const std::string& message) {
  std::lock_guard<std::mutex> lock(mutex_);

  if (level < level_) {
    return;
  }

  auto now = std::chrono::system_clock::now();
  auto now_c = std::chrono::system_clock::to_time_t(now);
  std::tm now_tm;

  localtime_r(&now_c, &now_tm);

  std::ostringstream timestamp_stream;
  timestamp_stream << std::put_time(&now_tm, "%Y-%m-%d %H:%M:%S");
  std::string timestamp = timestamp_stream.str();

  std::string str = level_to_string(level);
  std::string entry = timestamp + " [" + str + "] " + message;

  if (to_console_) {
    std::cout << entry << std::endl;
  }

  if (file_ && file_->is_open()) {
    *file_ << entry << std::endl;
  } else if (!path_.empty()) {
    if (to_console_) {
      std::cerr << timestamp << " [ERROR] Failed to write to log file: " << path_ << std::endl;
    } else {
      std::cerr << timestamp << " [ERROR] Failed to write to log file: " << path_ << std::endl;
      std::cerr << timestamp << " [ERROR] Original message: " << message << std::endl;
    }
  }
}

void Logger::debug(const std::string& message) {
  log(LogLevel::DEBUG, message);
}

void Logger::info(const std::string& message) {
  log(LogLevel::INFO, message);
}

void Logger::warn(const std::string& message) {
  log(LogLevel::WARNING, message);
}

void Logger::error(const std::string& message) {
  log(LogLevel::ERROR, message);
}
} // namespace libqm
