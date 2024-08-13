#include <iostream>
#include <fstream>   // std::ofstream
#include <vector>
#include <string>
#include <cassert>
#include <algorithm> // std::find
#include <chrono>    // std::strftime



#include "io.h"

rapidjson::Document get_json_from_file(std::string file_name) {
    std::ifstream ifs(file_name);
    if (!ifs) {
        std::cerr << "Error reading file ("
                  << file_name
                  << ")\n";
        exit(EXIT_FAILURE);
    }
    std::string file_contents(
        (std::istreambuf_iterator<char>(ifs)),
        std::istreambuf_iterator<char>()
    );
    ifs.close();

    rapidjson::Document rj_doc;
    rj_doc.Parse(file_contents.c_str());
    if (rj_doc.HasParseError()) {
        std::cerr << "File ("<< file_name << ")"
                  << "is not a valid JSON\n"
                  << "Error at offset "
                  << static_cast<unsigned>(rj_doc.GetErrorOffset())
                  << ": "
                  << rapidjson::GetParseError_En(rj_doc.GetParseError())
                  << "\n";
        exit(EXIT_FAILURE);
    }
    return rj_doc;
}

std::string get_time_stamp(){
    auto t = std::time(nullptr);
    auto tm = std::localtime(&t);
    char buff[16] = {0};
    std::strftime(buff, sizeof(buff), "%Y%m%d%H%M%S", tm);
    return std::string(buff);
}

InputParser::InputParser(int argc, char **argv) {
    for (int ii = 1; ii < argc; ii++){
        this->tokens.push_back(std::string(argv[ii]));
    }
}

const std::string& InputParser::getCmdOption(const std::string& option) const {
    std::vector<std::string>::const_iterator itr;
    itr = std::find(this->tokens.begin(), this->tokens.end(), option);

    if (itr != this->tokens.end() && ++itr != this->tokens.end()){
        return *itr;
    }

    static const std::string empty_string("");
    return empty_string;
}

bool InputParser::cmdOptionExists(const std::string& option) const {
    return std::find(this->tokens.begin(), this->tokens.end(), option)
           != this->tokens.end();
}




template <typename T>
void write_to_csv(
    std::string file_name,
    std::vector<std::vector<T>> vov,
    std::vector<std::string> v_names
    ) {
    assert(file_name != "");
    assert(!vov.empty());

    bool has_names = !v_names.empty();
    if (has_names) {
        assert(v_names.size() == vov.size());
    }

    std::ofstream ofs(file_name);
    std::cout << "writing to " << file_name << " ... ";
    if (vov.size() != 1) {
        for (size_t cc = 0; cc < vov.size(); cc++) {
            std::string cc_header = "v" + std::to_string(cc);
            if (has_names) {
                cc_header = v_names.at(cc);
            }
            ofs << cc_header << kCsv_sep;
        }
        ofs << "\n";
    }

    for (size_t rr = 0; rr < vov[0].size(); rr++) {
        for (size_t cc = 0; cc < vov.size(); cc++) {
            ofs << vov[cc][rr] << kCsv_sep;
        }
        ofs << "\n";
    }

    ofs.close();
    std::cout << " done.\n"; //<< file_name <<" closed.\n";
}
template void write_to_csv<int>(
    std::string file_name,
    std::vector<std::vector<int>> vov,
    std::vector<std::string> v_names = std::vector<std::string>{}
);
template void write_to_csv<float>(
    std::string file_name,
    std::vector<std::vector<float>> vov,
    std::vector<std::string> v_names = std::vector<std::string>{}
);
template void write_to_csv<double>(
    std::string file_name,
    std::vector<std::vector<double>> vov,
    std::vector<std::string> v_names = std::vector<std::string>{}
);
template void write_to_csv<std::string>(
    std::string file_name,
    std::vector<std::vector<std::string>> vov,
    std::vector<std::string> v_names = std::vector<std::string>{}
);

template <typename T>
void write_to_csv(
    std::string file_name,
    std::vector<T> v,
    std::string v_name
    ) {
    write_to_csv(
        file_name,
        std::vector<std::vector<T>>{v},
        std::vector<std::string>{v_name}
    );
}
template void write_to_csv<int>(std::string file_name, std::vector<int> v, std::string v_name = "");
template void write_to_csv<float>(std::string file_name, std::vector<float> v, std::string v_name = "");
template void write_to_csv<double>(std::string file_name, std::vector<double> v, std::string v_name = "");
template void write_to_csv<std::string>(std::string file_name, std::vector<std::string> v, std::string v_name = "");

template <typename T>
void print_v(std::vector<T> v, std::string v_name) {
    std::cout << "v(" << v_name << ")" << ":";
    for (T ee : v) {
        std::cout << ee << ", ";
    }
    std::cout << "\n";
}
template void print_v<int>(std::vector<int> v, std::string v_name = "");
template void print_v<float>(std::vector<float> v, std::string v_name = "");
template void print_v<double>(std::vector<double> v, std::string v_name = "");
template void print_v<std::string>(std::vector<std::string> v, std::string v_name = "");

template <typename T>
void print_vov( // print vector of vectors
    std::vector<std::vector<T>> vov,
    std::vector<std::string> v_names
    ) {
// void print_vov(auto vov) { // c++20
    bool has_names = !v_names.empty();
    if (has_names) {
        assert(v_names.size() == vov.size());
    }
    // std::cout << "v_names.size() : " << v_names.size() << "\n";
    int ii = 0;
    for (std::vector<T> vv : vov) {
        if (has_names) {
            print_v(vv, v_names[ii++]);
        } else {
            print_v(vv, std::to_string(ii++));
        }
    }
}
template void print_vov<int>(
    std::vector<std::vector<int>> vov,
    std::vector<std::string> v_names = std::vector<std::string>{}
    );
template void print_vov<float>(
    std::vector<std::vector<float>> vov,
    std::vector<std::string> v_names = std::vector<std::string>{}
    );
template void print_vov<double>(
    std::vector<std::vector<double>> vov,
    std::vector<std::string> v_names = std::vector<std::string>{}
    );
template void print_vov<std::string>(
    std::vector<std::vector<std::string>> vov,
    std::vector<std::string> v_names = std::vector<std::string>{}
    );