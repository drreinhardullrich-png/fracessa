#include <fracessa/logger.hpp>

Logger::Logger(const std::string& filename) {
    if (!filename.empty()) {
        _log_file.emplace(filename);
    }
}
