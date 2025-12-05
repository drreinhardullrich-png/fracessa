#ifndef CANDIDATE_H
#define CANDIDATE_H

//#include <vector>
//#include <string>

#include "Helper.hpp"
#include "Matrix.hpp"

class Candidate
{
    public:
        size_t candidate_id = 0;
        std::vector<rational> vector;
        uint64_t support;
        size_t support_size;
        uint64_t extended_support;
        size_t extended_support_size;
        size_t shift_reference;
        bool is_ess;
        ReasonEss reason_ess;
        rational payoff;
        double payoff_double;


        std::string to_string()
        {
            std::string str = "";
            str += std::to_string(candidate_id) + ";";
            for (auto x : vector)
                str += x.template convert_to<std::string>() + ",";
            str.pop_back();
            str += ";" + std::to_string(support) + ";";
            str += std::to_string(support_size) + ";";
            str += std::to_string(extended_support) + ";";
            str += std::to_string(extended_support_size) + ";";
            str += std::to_string(shift_reference) + ";";
            str += std::to_string(is_ess) + ";";
            str += std::to_string((int)reason_ess) + ";";
            str += payoff.template convert_to<std::string>() + ";";
            str += std::to_string(payoff_double);

            return str;
        }

        static std::string header()
        {
            return "candidate_id;vector;support;support_size;extended_support;extended_support_size;shift_reference;is_ess;reason_ess;payoff;payoff_double;";
        }
};

#endif // CANDIDATE_H
