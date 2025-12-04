#ifndef MATRIX_H
#define MATRIX_H

#include <vector>
#include <exception>
#include <sstream>

#include <fracessa/helper.hpp>

class matrix_exception: public std::exception
{
  virtual const char* what() const noexcept override
  {
  return "Exception while creating a matrix instance. Your input is not correct.";
  }
};


template <typename T>
class matrix
{
  static_assert(std::is_same<T, double>::value || std::is_same<T, rational>::value, "T is neither rational nor double!");

  private:

    std::vector<T> _matrix;
    size_t _rows;
    size_t _cols;
    size_t _length;
    bool _cs = false;
    //static matrixException _mex;

  public:

    bool is_cs() const {return _cs;}
    size_t rows() const {return _rows;}
    size_t cols() const {return _cols;}

    //*************************************constructors************************************************************
    matrix() {}

    matrix(size_t rows, size_t columns)
    {
      _rows = rows;
      _cols = columns;
      _length = rows * columns;
      _matrix.resize(_length);
    }

    //********************************Rule of Five************************************************************

    // Destructor
    ~matrix() = default;

    // Copy constructor
    matrix(const matrix&) = default;

    // Copy assignment operator
    matrix& operator=(const matrix&) = default;

    // Move constructor
    matrix(matrix&&) noexcept = default;

    // Move assignment operator
    matrix& operator=(matrix&&) noexcept = default;

    //********************************static construction methods*************************************************************************


    //full nxn-matrix given as vector or string of len n^2, circular symmetric optimizations NOT applied even if matrix is cs!

    template <class U>
    static inline matrix<rational> create_general(const std::vector<U> &elements)
    {
      size_t n;
      try {
        n = std::sqrt(elements.size());
      } catch (std::exception& e) {
        throw std::runtime_error("Number of elements not correct!");
      }

      matrix<rational> A = matrix<rational>(n, n);

      try {
        for (size_t i=0; i<elements.size(); i++) {
          A._matrix[i] = rational(elements[i]);
        }
      } catch (std::exception& e) {
        throw std::runtime_error("Could not create rational numbers from the elements!");
      }

      return A;
    }

    static inline matrix<rational> create_general(const std::string &elements)
    {
      std::vector<std::string> splitted;
      boost::split(splitted, elements, boost::is_any_of(","));
      return matrix<rational>::create_general(splitted);
    }


    //circular nxn-matrices given by the first row, i.e. a vector (or string) of len n has to be provided. checks if the matrix is cs automatically.
    //note that the first entry in the vector/string does not change the number or structure of solutions, it only changes the payoff.

    template <class U>
    static inline matrix<rational> create_circular(const std::vector<U> &first_row)
    {
      size_t n = first_row.size();

      matrix<rational> A = matrix<rational>(n, n);

      bool is_cs = true;
      for (size_t i=0; i<(n-1)/2; i++)
        if (first_row[i+1] != first_row[n-i-1]) {
          is_cs = false;
          break;
        }
      A._cs = is_cs;

      try {
        for (size_t i=0; i<n; i++)
          for (size_t j=0; j<n; j++)
            A(i, j) = rational(first_row[(j-i+n)%n]);
      } catch (std::exception& e) {
        throw std::runtime_error("Could not create rational number from first_row!");
      }

      return A;
    }

    static inline matrix<rational> create_circular(const std::string &first_row)
    {
      std::vector<std::string> splitted;
      boost::split(splitted, first_row, boost::is_any_of(","));
      return matrix<rational>::create_circular(splitted);
    }


    //circular nxn-matrices given by the first half-row starting at the second element, i.e. a vector (or string) of len floor(n/2) has to be provided. n has to be given bec. the vector/string alone is not sufficient to determine n.
    //as first entry in the entire row zero is prepended. this does not change the solutions, see above.
    //examples: for parameters 4, [7,13] the whole first row becomes [0, 7, 13, 7]
    //examples: for parameters 5, [7,13] the whole first row becomes [0, 7, 13, 13, 7]

    template <class U>
    static inline matrix<rational> create_circular_symmetric(size_t n, const std::vector<U> &half_row)
    {
      size_t half_row_size = half_row.size();

      if (half_row_size != n/2)
        throw std::runtime_error("Size of vector half_row is not feasible for given n!");

      std::vector<rational> first_row(n);

      first_row[0] = 0;
      try {
        if (n%2 == 0) {//even n
          first_row[n/2] = rational(half_row[half_row_size-1]);
          for (size_t i=0; i<n/2-1; i++) {
            first_row[i+1] = rational(half_row[i]);
            first_row[n-i-1] = rational(half_row[i]);
          }
        } else { //odd n
          for (size_t i=0; i<n/2; i++) {
            first_row[i+1] = rational(half_row[i]);
            first_row[n-i-1] = rational(half_row[i]);
          }
        }
      } catch (std::exception& e) {
        throw std::runtime_error("Could not create rational number from half_row!");
      }
      return matrix<rational>::create_circular(first_row);

    }

    static inline matrix<rational> create_circular_symmetric(size_t n, const std::string &half_row)
    {
      std::vector<std::string> splitted;
      boost::split(splitted, half_row, boost::is_any_of(","));
      return matrix<rational>::create_circular_symmetric(n, splitted);
    }


    //input format for cli interface. see documentation

    static inline matrix<rational> create_from_cli_string(const std::string &input)
    {
      std::vector<std::string> first_split;
      boost::split(first_split, input, boost::is_any_of("#"));
      if (first_split.size() != 2) {
        throw std::runtime_error("String for the matrix does not include '#' as a seperator between dimension and matrix, or several '#' were found!");
      }
      size_t n;
      try {
        n = std::stoull(first_split[0]);
      } catch (std::exception& e) {
        throw std::runtime_error("The given dimension could not be converted into an integer number!");
      }
      std::vector<std::string> second_split;
      boost::split(second_split, first_split[1], boost::is_any_of(","));

      if (second_split.size() == n/2) { //circular symmetric
        return matrix<rational>::create_circular_symmetric(n, second_split);
      } else if (second_split.size() == n*n) {
        return matrix<rational>::create_general(second_split);
      } else {
        throw std::runtime_error("The number of matrix-elements must either be floor(dimension/2) (for a cyclically symmetric matrix) or dimension^2 (for any other matrix)! Neither of that is the case!");
      }
    }

    //rest

    static inline matrix<T> zero_matrix(int rows, int cols)
    {
      matrix<T> result = matrix<T>(rows, cols);
      for (int i=0; i<rows; i++)
        for (int j=0; j<cols; j++)
          result(i,j) = 0;
      return result;
    }

    static inline matrix<T> identity_matrix(size_t n)
    {
      matrix<T> result = zero_matrix(n, n);
      for (size_t i=0; i<n; i++)
        result(i,i) = 1;
      return result;
    }

    //**********************************************access elements*****************************************************

    inline T& operator()(size_t row, size_t col)
    {
      return _matrix[row*_cols+col];
    }

    inline const T& operator()(size_t row, size_t col) const
    {
      return _matrix[row*_cols+col];
    }


    //***************************************************everything else****************************************************

    inline std::string to_string() const
    {
      std::ostringstream oss;
      for (size_t i=0;i<_rows;i++) {
        for (size_t j=0; j<_cols;j++)
          oss << (*this)(i,j) << ",";
        oss << std::endl;
      }
      return oss.str();
    }

    inline matrix<double> to_double() const
    {
      static_assert(std::is_same<T, rational>::value, "T is not rational!");

      matrix<double> A = matrix<double>(_rows, _cols);
      for (size_t i=0; i<_rows; i++)
        for (size_t j=0; j<_cols;j++) {
          A(i,j) = static_cast<double>(_matrix[i*_cols+j]);
        }
      return A;
    }

    inline void get_le_matrix(uint64_t support, size_t support_size, matrix<T> &le_matrix)
    {
      size_t row = 0;
      size_t column = 0;

      for (size_t i=0; i<_rows; i++)
        if ((support & (1ull << i)) != 0)
        {
          column = 0;
          for (size_t j=0; j<_cols; j++)
            if ((support & (1ull << j)) != 0)
            {
              le_matrix(row, column) = _matrix[i*_cols+j];
              column++;
            }
          le_matrix(row, column) = -1;
          le_matrix(row, column+1) = 0;
          row++;

        }
      for (size_t i=0; i<support_size; i++)
        le_matrix(support_size, i) = 1;

      le_matrix(support_size, support_size) = 0;
      le_matrix(support_size, support_size+1) = 1;
    }

    inline bool is_positive_definite()
    {
      static_assert(std::is_same<T, rational>::value, "T is not rational!");
      //ldlt-decomposition for rational numbers
      //      //from sympy/matrices/dense.py
      //      D = zeros(self.rows, self.rows)
      //          L = eye(self.rows)
      //          for i in range(self.rows):
      //              for j in range(i):
      //                  L[i, j] = (1 / D[j, j])*(self[i, j] - sum(L[i, k]*L[j, k]*D[k, k] for k in range(j)))
      //              D[i, i] = self[i, i] - sum(L[i, k]**2*D[k, k] for k in range(i))
      //              return self._new(L), self._new(D)

      int n = _rows;
      matrix<rational> D = zero_matrix(n, n);
      matrix<rational> L = identity_matrix (n);

      for (int i = 0; i < n; i++) {
        for (int j = 0; j < i; j++) {
          rational aSum = rational(0);
          for (int k = 0; k < j; k++)
            aSum += L(i,k) * L(j,k) * D(k,k);
          L(i,j) = (1 / D(j,j)) * ((*this)(i,j) - aSum);
        }
        rational bSum = rational(0);
        for (int k = 0; k < i; k++)
          bSum += L(i,k) * L(i,k) * D(k,k);
        D(i,i) = (*this)(i,i) - bSum;
        if (D(i,i)<=rational(0))
          return false;
      }
      return true;
    }

    inline bool is_positive_definite_double()
    {
      static_assert(std::is_same<T, double>::value, "T is not double!");
      // cholesky-decomposition for double!
      //      //from sympy/matrices/dense.py
      //      L = zeros(self.rows, self.rows)
      //          for i in range(self.rows):
      //              for j in range(i):
      //                  L[i, j] = (1 / L[j, j])*(self[i, j] -sum(L[i, k]*L[j, k] for k in range(j)))
      //              L[i, i] = sqrt(self[i, i] -sum(L[i, k]**2 for k in range(i)))
      //          return self._new(L)

      int n = _rows;
      matrix<double> L = zero_matrix(n, n);

      for (int i = 0; i < n; i++) {
        for (int j = 0; j < i; j++) {
          double aSum = 0.;
          for (int k = 0; k < j; k++)
            aSum += L(i,k) * L(j,k);
          L(i,j) = (1 / L(j,j)) * ((*this)(i,j) - aSum);
        }
        double bSum = 0.;
        for (int k = 0; k < i; k++)
          bSum += L(i,k) * L(i,k);
        double x = (*this)(i,i) - bSum;
        if (x< 0.01)
          return false;
        L(i,i) = sqrt(x);
      }
      return true;
    }

    inline bool is_copositive()
    {
      static_assert(std::is_same<T, rational>::value, "T is not rational!");

      for (uint64_t support = 1ull; support < (1ull << _rows); support++) { //iterate all subsets without the empty set
        size_t n = support_size_from_int(support);
        matrix<T> A = matrix<T>(n,n);

        size_t row = 0;
        size_t column = 0;
        for (size_t i=0; i<_rows; i++) {
          if ((support & (1ull << i)) != 0) {
            column = 0;
            for (size_t j=0; j<_cols; j++)
              if ((support & (1ull << j)) != 0) {
                A(row,column) = (*this)(i,j);
                column++;
              }
            row++;
          }
        }
        if (A.determinant() <= 0 && A.adjugate().greater_zero()) {
          return 0;
        }

      }
      return 1;
    }

    inline matrix<T> clone()
    {
      matrix<T> A = matrix<T> (_rows, _cols);
      for (size_t i=0; i<_rows; i++)
        for (size_t j=0; j<_cols;j++) {
          A(i,j) = _matrix[i*_cols+j];
        }
      return A;
    }

    inline T determinant()
    {
      static_assert(std::is_same<T, rational>::value, "T is not rational!");
      //lu-decomposition, crout-algorithm from http://algorithm.narod.ru/ln/crout.c, only for rational bec. check for zero!

      matrix<T> a = this->clone();
      size_t n = _rows;
      std::vector<size_t> indx = std::vector<size_t>(n);
      std::vector<T> vv = std::vector<T>(n);
      int d = 1;

      size_t i, imax, j, k;
      T big, sum, temp;

      /* search for the largest element in each row; save the scaling in the
      temporary array vv and return zero if the matrix is singular */
      for(i=0; i<n; i++) {
        big=0.;
        for(j=0;j<n;j++) if((temp=abs(a(i,j)))>big) big=temp;
        if(big==0) return(0);
        vv[i]=big;
       }
       /* the main loop for the Crout's algorithm */
       for(j=0;j<n;j++) {
        /* this is the part a) of the algorithm except for i==j */
        for(i=0;i<j;i++) {
          sum=a(i,j);
          for(k=0;k<i;k++) sum-=a(i,k)*a(k,j);
          a(i,j)=sum;
        }
        /* initialize for the search for the largest pivot element */
        big=0;imax=j;
        /* this is the part a) for i==j and part b) for i>j + pivot search */
        for(i=j;i<n;i++) {
          sum=a(i,j);
          for(k=0;k<j;k++) sum-=a(i,k)*a(k,j);
          a(i,j)=sum;
          /* is the figure of merit for the pivot better than the best so far? */
          if((temp=vv[i]*abs(sum))>=big) {big=temp;imax=i;}
        }
        /* interchange rows, if needed, change parity and the scale factor */
        if(imax!=j) {
          for(k=0;k<n;k++) {temp=a(imax,k);a(imax,k)=a(j,k);a(j,k)=temp;}
          d=-d;vv[imax]=vv[j];
        }
        /* store the index */
        indx[j]=imax;
        if(a(j,j)==0) return(0);
        /* finally, divide by the pivot element */
        if(j<n-1) {
          temp=1/a(j,j);
          for(i=j+1;i<n;i++) a(i,j)*=temp;
        }
      }
      T res = d;
      for(j=0; j<n; j++) res *= a(j,j);
      return(res);
    }

    inline matrix<T> transpose()
    {
      matrix<T> A = matrix<T>(_rows,_cols);
      for (size_t i=0; i<_rows; i++)
        for (size_t j=0; j<_cols; j++)
          A(i,j) = (*this)(j,i);
      return A;
    }

    inline matrix<T> develop_by(size_t row, size_t col)
    {
      matrix<T> A=matrix<T>(_rows-1,_cols-1);
      for (size_t i=0; i<_rows-1; i++)
        for (size_t j=0; j<_cols-1; j++) {
          size_t row_index_real = (i>=row) ? i+1 : i;
          size_t col_index_real = (j>=col) ? j+1 : j;
          A(i,j) = (*this)(row_index_real,col_index_real);
        }
      return A;
    }

    inline matrix<T> adjugate()
    {
      matrix<T> A = matrix<T>(_rows,_cols);

      if (_rows == 1) {
        A(0,0) = 1;
        return A;
      }

      for (size_t i=0; i<_rows; i++) {
        for (size_t j=0; j<_cols; j++) {
          A(i,j) = std::pow(-1, i+j) * (this->develop_by(i,j)).determinant();
        }

      }
      return A.transpose();
    }

    inline bool greater_zero()
    {
      for (size_t i=0; i<_rows; i++)
        for (size_t j=0; j<_cols; j++)
          if ((*this)(i,j) <= 0)
            return false;
      return true;
    }

};

#endif // MATRIX_H

