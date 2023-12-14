//
// Created by mafu on 12/13/2023.
//

#pragma once

#include <fstream>
#include <iostream>
#include <map>
#include <string>

class config_parser {
    std::map<std::string, std::string> config_data_;

    static std::string trim(const std::string &s) {
        auto start = s.begin();
        while (start != s.end() && std::isspace(*start)) {
            ++start;
        }

        auto end = s.end();
        do {
            --end;
        } while (end != start && std::isspace(*end));

        return {start, end + 1};
    }

public:
    explicit config_parser(const std::string &filename) {
        std::ifstream file(filename);
        if (!file.is_open()) {
            std::cerr << "Error: Cannot open config file: " << filename << std::endl;
            throw std::runtime_error("Cannot open config file");
        }

        std::string line;
        while (std::getline(file, line)) {
            line = line.substr(0, line.find(';'));
            const size_t pos = line.find('=');
            if (pos != std::string::npos) {
                std::string key = line.substr(0, pos);
                std::string value = line.substr(pos + 1);
                config_data_[trim(key)] = trim(value);
            }
        }

        file.close();
    }

    [[nodiscard]] std::string get(const std::string &key) const {
        try {
            return config_data_.at(key);
        }
        catch (const std::out_of_range &) {
            return "";
        }
    }

    [[nodiscard]] int get_int(const std::string &key) const {
        try {
            return std::stoi(config_data_.at(key));
        }
        catch (const std::out_of_range &) {
            return 0;
        }
    }

    [[nodiscard]] double get_double(const std::string &key) const {
        try {
            return std::stod(config_data_.at(key));
        }
        catch (const std::out_of_range &) {
            return 0.0;
        }
    }

    [[nodiscard]] bool get_bool(const std::string &key) const {
        try {
            const std::string value = config_data_.at(key);
            return value == "True" || value == "true" || value == "1";
        }
        catch (const std::out_of_range &) {
            return false;
        }
    }
};

