#ifndef LOGGER_H
#define LOGGER_H

#include <fstream>
#include <string>
#include <optional>
#include <memory>

class Logger {
public:
    Logger(const std::string& filename);

    template<typename T>
    Logger& operator<<(const T& value) {
        if (_log_file) {
            *_log_file << value;
        }
        return *this;
    }

    // Overload for std::endl and other manipulators
    Logger& operator<<(std::ostream& (*pf)(std::ostream&)) {
        if (_log_file) {
            *_log_file << pf;
        }
        return *this;
    }

    bool is_enabled() const {
        return _log_file.has_value();
    }

private:
    std::optional<std::ofstream> _log_file;
};

#endif // LOGGER_H
