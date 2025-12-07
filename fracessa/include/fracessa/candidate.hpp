#ifndef CANDIDATE_H
#define CANDIDATE_H

#include <fracessa/matrix.hpp>
#include <fracessa/bitset64.hpp>

class candidate
{
    public:
        size_t candidate_id = 0;
        RationalVector vector;
        bitset64 support;
        size_t support_size;
        bitset64 extended_support;
        size_t extended_support_size;
        size_t shift_reference;
        bool is_ess;
        std::string stability;
        rational payoff;
        double payoff_double;


        std::string to_string()
        {
            std::string str = "";
            str += std::to_string(candidate_id) + ";";
            for (Eigen::Index i = 0; i < vector.size(); i++) {
                str += vector(i).template convert_to<std::string>() + ",";
            }
            if (vector.size() > 0)
                str.pop_back();
            str += ";" + support.to_string() + ";";
            str += std::to_string(support_size) + ";";
            str += extended_support.to_string() + ";";
            str += std::to_string(extended_support_size) + ";";
            str += std::to_string(shift_reference) + ";";
            str += std::to_string(is_ess) + ";";
            str += stability + ";";
            str += payoff.template convert_to<std::string>() + ";";
            str += std::to_string(payoff_double);

            return str;
        }

        static std::string header()
        {
            return "candidate_id;vector;support;support_size;extended_support;extended_support_size;shift_reference;is_ess;stability;payoff;payoff_double;";
        }
};

#endif // CANDIDATE_H
