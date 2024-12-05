// writeresults.h

#ifndef WRITERESULTS_H
#define WRITERESULTS_H


class writeresults {
public:
    //template <typename... Args>
    //void save_bin_data(const char* , const Args&... );

    void save_bin_data(const char* , const char** , const double* , int);
//private:
//    template <typename T>
//    void writeHeader(std::ostream& file, const T& arg);

//    template <typename T, typename... Args>
//    void writeHeader(std::ostream& file, const T& arg, const Args&... args);
};

#endif // WRITERESULTS_H

