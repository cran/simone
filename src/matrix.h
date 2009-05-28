//  matrix.h
// 
//  Copyright (C) 2008 Laboratoire Statistique & Genome
//
//  This program is free software; you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation; either version 2 of the License, or (at
//  your option) any later version.
// 
//  This program is distributed in the hope that it will be useful, but
//  WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//  General Public License for more details.
// 
//  You should have received a copy of the GNU General Public License
//  along with this program; if not, write to the Free Software
//  Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
// 
//

//! \file    matrix.h
//! \author  Gilles Grasseau
//! \brief   Basic tools for matrix
//! \version 0.1
//! \date    June 2008

# include <iostream>
# include <cstring>
# include <cmath>

# include "util.h"

# ifndef _MATRIX_H_
# define _MATRIX_H_

static const char *not_same_dims   = "Matrix don't have the same dimensions";
static const char *not_dims_match  = "Matrix dimensions don't match";
static const char *not_scalar      = "Matrix can't be trasformed to scalar";
static const char *not_implemented = "Operation/function not implemented";

//! \cond
#define EQUAL_DIM( X ) ( (row_ == (X).row_) && (col_ == (X).col_) ) 

#define EQUAL_DIM_MAT( X, Y ) ( ((X).getRow() == (Y).getRow()) \
                             && ((X).getCol() == (X).getCol()) )

#define EQUAL_COL_ROW( X ) ( (col_ == (X).row_) )  
//! \endcond

//! \class Matrix matrix.h
//! \brief Class used to make matrix computations
//!  The matrix elements are row-major ordered i.e. A(i,k) and A(i+1,k) are 
//!  contiguous in the memory. 
class Matrix {

  private :

  size_t row_; 
  size_t col_;

  struct s_data {
    // Reference number on matrix array matrix_
    int ref_;
    // Row ordered first
    double *matrix_;
  } *link_;
 
  struct s_data *newLink(size_t size_ ){
    LOGInFct( std::cout );
    struct s_data *link = new struct s_data;
    link -> matrix_ = new double [size_];
    link -> ref_    = 1;
    LOGMSG( 4, std::cout, "New link", (long int) link->matrix_ );

    LOGOutFct( std::cout );
    return(link);
  }

  struct s_data *newLink( double *ptr_ ){
    LOGInFct( std::cout );
    struct s_data *link = new struct s_data;
    link -> matrix_ = ptr_;
    link -> ref_    = 1;
    LOGMSG( 4, std::cout, "New link", (long int) link->matrix_ );
    LOGOutFct( std::cout );
    return(link);
  }

  void deleteLink( ){
    LOGInFct( std::cout );
    if ( (-- (link_-> ref_)) == 0 ) {
      LOGMSG( 4, std::cout, "Freeing link", (long int) link_->matrix_ );
      delete [] (link_->matrix_);
      delete link_;
      link_ = 0;
    }
    LOGOutFct( std::cout );
  }

  void fillData( double a ) {     
    for (size_t i=0; i < row_*col_; i++) {
      link_->matrix_[i] = a;
    }
  } 

  void fillData( const double *ptr ) {   

    memcpy( link_->matrix_, ptr, row_*col_*sizeof(double) );
  } 

  public :

    //! \brief Constructor - The matrix is not created.
    Matrix() {
      row_ = 0; col_ = 0;
      link_ = 0;
    }

    //! \brief Constructor - The matrix A(row,col) contains indeterminated
    //!  values.
    //! \param row Number of rows.
    //! \param col Number of lines.
    Matrix( const size_t row, const size_t col) {
      row_ = row; col_ = col;
      link_ = newLink( row_*col_ );
    }

    //! \brief Constructor - The matrix  A(row,col) is filled with a value.
    //! \param row Number of rows.
    //! \param col Number of lines.
    //! \param a   Value to fill the matrix with.
    Matrix( const size_t row, const size_t col, double a) {
      row_ = row; col_ = col;
      link_ = newLink( row_*col_ );
      fillData( a );
    }

    //! \brief Constructor -  The matrix  A(row,col) is filled with
    //!        the "ptr" pointer contain.
    //! \param row  Number of rows.
    //! \param col  Number of lines.
    //! \param ptr  Floating point array.
    //! \param copy Copy the  Floating point "ptr" array if true, 
    //!             else "ptr" is directly assign to the internal matrix data.
    //!             In "false" case do not use the pointer for other 
    //!              computations  
    Matrix( const size_t row, const size_t col, double *ptr, 
	    const bool copy=true ) {
      LOGInFct( std::cout );
      row_ = row; col_ = col;
      if ( copy ) {
	link_ = newLink(row_*col_) ;
	fillData( ptr );
      } else {
	link_ = newLink( ptr );
      }
      LOGOutFct( std::cout );
    }


    //! Copy constructor
    //! \param A The new matrix shared the same data with matrix A.
    //!        The new matrix is an alias of A.
    //!        Do not used this constuctor to duplicate A. 
    //! \code
    //!  B = Matrix(A);        // Forbiden to copy matrix
    //!  Matrix C = A;         // Forbiden to copy matrix
    //! \endcode
    //! See getCopy() to duplicate Matrix objects.
    Matrix( const Matrix& A ) {
      LOGInFct( std::cout );
      row_ = A.row_;
      col_ = A.col_;
      A.link_-> ref_++;
      link_ = A.link_;
      LOGOutFct( std::cout );
    }

    //! Destructor
    ~Matrix() {
      LOGInFct( std::cout );
      deleteLink();
      LOGOutFct( std::cout );
    }

    //! Assignement operator
    Matrix& operator = ( const Matrix& A) {
      LOGInFct( std::cout );
      row_ = A.row_;
      col_ = A.col_;
      deleteLink();

      link_ = A.link_;
      link_ -> ref_++;
      LOGOutFct( std::cout );

      return(*this);
    }

    //! \cond  Symetric operator
    friend Matrix inline operator+(const double a, const Matrix& A );
    friend Matrix inline operator-(const double a, const Matrix& A );
    friend Matrix inline operator*(const double a, const Matrix& A );
    //! \endcond

    //! \brief Compute the maximum between the A matrix elements and the scalar 
    //! value a.
    //! \param A Matrix.
    //! \param a Scalar.
    //! \return Result Matrix.
    friend inline Matrix max( const Matrix& A, const double  a  );
    //! \cond  Symetric operator
    friend inline Matrix max( const double  a, const Matrix& A );
    //! \endcond

    //! \brief Compute the maximum element-by-element between the A matrix 
    //! and the B matrix elements. 
    //! \param A Matrix.
    //! \param B Matrix.
    //! \return Result Matrix.
    friend inline Matrix max( const Matrix& A, const Matrix& B );

    //! \brief Compute the minimum between the A matrix elements and the scalar 
    //! value a.
    friend inline Matrix min( const Matrix& A, const double  a );
    //! \cond  Symetric operator
    friend inline Matrix min( const double  a, const Matrix& A );
    //! \endcond

    //! \brief Compute the minimun element-by-element between the A matrix 
    //! and the B matrix elements. 
    //! \param A Matrix.
    //! \param B Matrix.
    //! \return Result Matrix.
    friend inline Matrix min( const Matrix& A, const Matrix& B );

    //! \brief Compute the sum of all A matrix elements. 
    //!        sum = sum( A(i,j), 0 < i < row, 0 < j < col) 
    //! \param B Matrix.
    //! \return The sum of all elements.
    friend inline double sum( const Matrix& A );

    //! \brief Compute the log() of A matrix elements. 
    //!        X(i,j) = log( A(i,j) ); 0 <= i < row, 0 <= j < col 
    //! \param A Matrix.
    //! \return log(A) Matrix.
    friend inline Matrix log( const Matrix& A );

    // Operators 

    //! \brief Type cast operator from a matrix (1,1) to a scalar.
    inline operator double(void) {
      double a=0;
      LOGInFct( std::cout );
      if ((row_ == 1) && (col_ == 1)) {
	a = link_->matrix_[0];
      } else 
	LOGMSG( 1, std::cout, not_scalar, "" );
      LOGOutFct( std::cout );
      return( a );
    }

    //! \brief Indexing of the matrix elements
    //! \param i Row indice.
    //! \param j Column indice.
    //! \return Return A(i,j) value.
    double inline &operator()(size_t i, size_t j) { 
      return ( link_->matrix_[j*row_+i] ); 
    };

    //! \brief Linear indexing of the matrix array
    //! Remark: The matrix elements are row-major ordered.
    //! \param i Indice, 0 <= i < row*col .
    //! \return Return A(i) scalar.
    double inline &operator()(size_t i) { 
      return ( link_->matrix_[i] ); 
    };

    //!cond
    //! Functionnal operators for "const Matrix"
    double inline operator()(size_t i, size_t j) const { // For "const Matrix" 
      return ( link_->matrix_[j*row_+i] ); 
    };    
    double inline operator()(size_t i) const { // For "const Matrix" 
      return ( link_->matrix_[i] ); 
    };    
    // \endcond
    
    // \cond
    //! Unitary operator
    Matrix inline operator+( );
    Matrix inline operator-( );
    // \endcond

    //! \brief Addition operator of a Matrix with a scalar.
    //! All Matrix elements are added with "a" value.
    //! \param a Scalar to add.
    //! \return Matrix 
    Matrix inline operator+(const double a );
    //! \brief Subtraction operator of a Matrix with a scalar.
    //! All Matrix elements are substracted with "a" value.
    //! \param a Scalar to substract.
    //! \return Matrix 
    Matrix inline operator-(const double a );
    //! \brief Multiply operator of a Matrix with a scalar.
    //! All Matrix elements are multiplied with "a" value.
    //! \param a Scalar to multiply.
    //! \return Result Matrix 
    Matrix inline operator*(const double a );
    //! \brief Divide operator of a Matrix with a scalar.
    //! All Matrix elements are divided with "a" value.
    //! \param a Scalar to divide.
    //! \return Result Matrix 
    Matrix inline operator/(const double a );

    //! \brief Add element-by-element the two matrices.
    //! \param A Matrix to add.
    //! \return Result Matrix  
    Matrix inline operator+(const Matrix& A);
    //! \brief Substract element-by-element the two matrices.
    //! \param A Matrix to substract.
    //! \return Result Matrix  
    Matrix inline operator-(const Matrix& A);
    //! \brief Mutiply element-by-element the two matrices.
    //! \param A Matrix to mutiply.
    //! \return Result Matrix  
    Matrix inline operator|(const Matrix& A);
    //! \brief Matrices multiplication.
    //! \param A Matrix to multiply.
    //! \return Result Matrix  
    Matrix inline operator*(const Matrix& A);

    //! \brief Duplicate the Matrix.
    //! Perform a deep copy of the object
    //! \return Return a copy of the object.  
    Matrix getCopy();

    //! \brief Get the number of rows in Matrix.
    //! \return Return the row number.  
    size_t inline getRow( ) const { return ( row_ ); }; 
    //! \brief Get the number of rows in Matrix.
    //! \return Return the row number.  
    size_t inline getCol( ) const { return ( col_ ); }; 

    //! \brief Get the matrix array.
    //! The content of the pointer cannot be modified. 
    //! \return Return a pointer to the data array.  
    const double *getMatrix() const { return( link_->matrix_ ) ; }

    //! \brief Tranpose the matrix. 
    //! \param A Matrix to multiply.
    //! \return Return a new Matrix transposed.  
    Matrix inline t(void);

    //! \brief Get the row vector i.
    //! \param i Row index.
    //! \return Return the vector as Matrix object.  
    Matrix inline getRowVector( const size_t i ); 
    //! \brief Get the column vector j.
    //! The content of the pointer cannot be modified. 
    //! \param i Column index.
    //! \return Return the vector as Matrix object.  
    Matrix inline getColVector( const size_t j ); 
    //! \brief Get the upper triangular matrix.
    //! Return a Matrix (row, col):
    //! - if (diag = false):
    //!   - A(i,j), 0 <= i < row, i < j < col
    //!   - zero for other matrix elements.
    //! - else (diag = true):
    //!   - A(i,j), 0 <= i < row, i <= j < col
    //!   - zero for other matrix elements.
    //! \param diag Specify if the diagonal is included
    //!        (default diag=false).
    //! \return Return the vector as Matrix object.  
    Matrix inline getUpperTriang( const bool diag=false );
    //! \brief Get the lower triangular matrix.
    //! Return a Matrix (row, col)
    //! if (diag = false):
    //!  - A(i,j), 0 <= i < row, 0 <= j < i
    //!  - zero for other matrix elements.
    //! else (diag = true):
    //!  - A(i,j), 0 <= i < row, 0 <= j <= i
    //!  - zero for other matrix elements.
    //! \param diag Specify if the diagonal is included
    //!        (default diag=false).
    //! \return Return the vector as Matrix object.  
    Matrix inline getLowerTriang( const bool diag=false );

    //! \brief Conditionnal assignment
    //! \param test Logical test to evaluate.
    //! \param a    Value to test with matrix elements.
    //! \param tvalue If the test is true set the corresponding
    //!               matrice element to "tvalue".
    //! \param fvalue If the test is false set the corresponding
    //!               matrice element to "fvalue".
    //! \return Return the Matrix result.
    Matrix inline Compare( const char   *str_, const double a, 
		    const double t_val, const double f_val );
    //! \brief Conditionnal assignment
    //! \param test Logical test to evaluate. Only "==" available.
    //! \param B    Matrix. Evaluate the test element-by-element.
    //! \param tvalue If the test is true set the corresponding
    //!               matrice element to "tvalue".
    //! \param fvalue If the test is false set the corresponding
    //!               matrice element to "fvalue".
    //! \return Return the Matrix result.
    Matrix inline Compare( const char   *str_, const Matrix& A, 
		    const double t_val, const double f_val );

    //! \brief Display the matrix.
    //! \param name Matrix name.
    void print( const char *name );
};


//   -------------------------------------
//
//   Define inlined functions and operators
//
//   -------------------------------------

//
//   Operators
//

Matrix inline Matrix::operator+( ) {

  double *tab = new double [ row_ * col_ ];
  double *mat_ = link_->matrix_;

  for (size_t i=0; i < row_*col_; i++) {
    tab[i] =  mat_[i];
  }
  Matrix tmp(row_, col_, tab, false);

  return ( tmp );
}

Matrix inline Matrix::operator-( ) {

  double *tab = new double [ row_ * col_ ];
  double *mat_ = link_->matrix_;

  for (size_t i=0; i < row_*col_; i++) {
    tab[i] = - mat_[i];
  }
  Matrix tmp(row_, col_, tab, false);

  return ( tmp );
}

Matrix Matrix::operator+( const double a ) {

  double *tab = new double [ row_ * col_ ];
  double *mat_ = link_->matrix_;

  for (size_t i=0; i < row_*col_; i++) {
    tab[i] = mat_[i]+a;
  }
  Matrix tmp(row_, col_, tab, false);

  return ( tmp );
}

Matrix  Matrix::operator+(const Matrix& A) {
    
  Matrix tmp( row_,col_);

  double *tmp_tab = tmp.link_->matrix_; 
  double *_tab    = link_->matrix_;
  double *A_tab   = A.link_->matrix_;
 
  if( EQUAL_DIM( A ) ) {
    
    for (size_t i=0; i < row_*col_; i++) {
      tmp_tab[i] = _tab[i] +  A_tab[i];
    }
  }
  else {
    LOGMSG( 1, std::cout, not_same_dims, "" );
  }
  return ( tmp );

}

Matrix Matrix::operator-( const double a ) {

  double *tab = new double [ row_ * col_ ];
  double *mat_ = link_->matrix_;

  for (size_t i=0; i < row_*col_; i++) {
    tab[i] = mat_[i] - a;
  }
  Matrix tmp(row_, col_, tab, false);

  return ( tmp );
}

Matrix inline Matrix::operator-(const Matrix& A) {
    
  Matrix tmp( row_,col_);

  double *tmp_tab = tmp.link_->matrix_; 
  double *_tab    = link_->matrix_;
  double *A_tab   = A.link_->matrix_;

  if( EQUAL_DIM( A ) ) {
    for (size_t i=0; i < row_*col_; i++) {
      tmp_tab[i] = _tab[i] -  A_tab[i];
    }
  }
  else {
    LOGMSG( 1, std::cout, not_same_dims, "" );
  }
  return ( tmp );

}

Matrix Matrix::operator*( double a ) {

  double *tab = new double [ row_ * col_ ];
  double *mat_ = link_->matrix_;

  for (size_t i=0; i < row_*col_; i++) {
    tab[i] = mat_[i] * a;
  }
  Matrix tmp(row_, col_, tab, false );

  return ( tmp );
}

Matrix inline Matrix::operator*(const Matrix& A) {
    
  Matrix tmp(row_, A.col_) ; 
  double *tmp_tab = tmp.link_->matrix_; 
  double *_tab    = link_->matrix_;
  double *A_tab   = A.link_->matrix_;
 
  LOGInFct( std::cout );

  if( EQUAL_COL_ROW( A ) ) {

    size_t A_row = A.row_;

#ifndef OPT_MULT
    // Initial loops
    double sum = 0.0;
    for (size_t i=0; i < row_; i++) {
      for (size_t j=0; j < A.col_; j++) {
	sum = 0.0;
	for (size_t k=0; k < col_; k++) {
	  sum += _tab[k*row_+i] * A_tab[j*A_row + k];
	}
	tmp_tab[j*row_+i] = sum;
      }
    }
#else
    double *sum_tab = new double[row_]; 
    // Loops i and k have switched for performance reason.
    for (size_t j=0; j < A.col_; j++) {
      memset(sum_tab, 0, sizeof(double) * row_ );
      for (size_t k=0; k < col_; k++) {
	for (size_t i=0; i < row_; i++) {
	  sum_tab[i] += _tab[k*row_+i] * A_tab[j*A_row + k];
	}
      }
      for (size_t i=0; i < row_; i++) {
	tmp_tab[j*row_+i] = sum_tab[i];
      }
    }
    free( sum_tab );
#endif
  }
  else {
    LOGMSG( 1, std::cout, not_dims_match, "" );
  }
  LOGOutFct( std::cout );
  return ( tmp );
}

Matrix inline Matrix::operator/( double a ) {

  double *tab = new double [ row_ * col_ ];
  double *mat_ = link_->matrix_;

  for (size_t i=0; i < row_*col_; i++) {
    tab[i] = mat_[i]/a;
  }
  Matrix tmp(row_, col_, tab, false);

  return ( tmp );
}

Matrix inline Matrix::operator|(const Matrix& A) {
    
  Matrix tmp(row_, col_) ; 

  double *tmp_tab = tmp.link_->matrix_; 
  double *_tab    = link_->matrix_;
  double *A_tab   = A.link_->matrix_;

  LOGInFct( std::cout );
  if( EQUAL_DIM( A ) ) {
    for (size_t i=0; i < row_*col_; i++) {
      tmp_tab[i] = _tab[i] * A_tab[i]; 
    }
  }
  else {
    LOGMSG( 1, std::cout, not_dims_match, "" );
  }
  LOGOutFct( std::cout );

  return ( tmp );
}

//
//    Member functions
//    ----------------

Matrix inline Matrix::t() {
    
  size_t t_row = col_;
  size_t t_col = row_;
  Matrix tmp(t_row, t_col) ; 
  double *tmp_tab = tmp.link_->matrix_; 
  double *_tab    = link_->matrix_;

  LOGInFct( std::cout );

  for (size_t j=0; j < t_col; j++) {
    for (size_t i=0; i < t_row; i++) {
      tmp_tab[j*t_row+i] = _tab[i*row_+j];
    }
  }

  LOGOutFct( std::cout );

  return ( tmp );
}

Matrix inline Matrix::getRowVector( const size_t i ) {

  const double* mat_ = getMatrix();
  double *tmp_mat    = new double [col_];
  for (size_t j=0; j < col_; j++) {
    tmp_mat[j] = mat_[j*row_+i];
  }
  
  Matrix tmp( 1, col_, tmp_mat, false );
  return( tmp );
}

Matrix inline Matrix::getColVector( const size_t j ) {

  const double* mat_ = getMatrix();
  double *tmp_mat    = new double [row_];

  memcpy ( tmp_mat, &mat_[j*row_], row_*sizeof(double));
  
  Matrix tmp( row_, 1, tmp_mat, false );
  return( tmp );
}

Matrix inline Matrix::getUpperTriang(const bool diag) {
  Matrix tmp = getCopy();
  double *mat_= tmp.link_->matrix_;
  size_t shift;

  if( diag ) 
    shift = 1;
  else
    shift = 0;

  for (size_t j=0; j < col_; j++) {
    for (size_t i=j+shift; i < row_; i++) {
      mat_[j*row_+i] = 0;
    }
  }
  return( tmp );
}

Matrix inline Matrix::getLowerTriang(const bool diag) {
  Matrix tmp = getCopy();
  double *mat_= tmp.link_->matrix_;
  size_t shift;
  if( diag ) 
    shift = 1;
  else
    shift = 0;
  for (size_t j=shift; j < col_; j++) {
    for (size_t i=0; i < j-1; i++) {
      mat_[j*row_+i] = 0;
    } 
  } 
  return( tmp );
}

Matrix inline Matrix::Compare( const char *str_,   const Matrix& A, 
			const double t_val, const double f_val ) {

  std::string str(str_);
  const double *mat_  = getMatrix();
  const double *mat_A = A.getMatrix();
  double *tmp_mat = new double[row_*col_];

  if( EQUAL_DIM( A ) ) {
    if ( str.compare("==") == 0){
      for (size_t i=0; i < row_*col_; i++) {
	if ( mat_[i] ==  mat_A[i]) 
	  tmp_mat[i] = t_val;
	else
	  tmp_mat[i] = f_val;
      }
    }
    else {
      LOGMSG( 1, std::cout, not_implemented, "" );
    }
  } else {
    LOGMSG( 1, std::cout, not_same_dims, "" );
  }
  Matrix tmp( col_, row_, tmp_mat,  false);
  return ( tmp );

}

Matrix inline Matrix::Compare( const char *str_,   const double test_val, 
			const double t_val, const double f_val ) {

  std::string str(str_);
  const double *mat_  = getMatrix();
  double *tmp_mat = new double[row_*col_];

  if ( str.compare("==") == 0){
    for (size_t i=0; i < row_*col_; i++) {
      if ( mat_[i] == test_val ) 
	tmp_mat[i] = t_val;
      else
	tmp_mat[i] = f_val;
    }
  }
  else {
    LOGMSG( 1, std::cout, not_implemented, "" );
  }
  Matrix tmp( col_, row_, tmp_mat, false);
  return ( tmp );

}

//
//   Friend functions and operators
// 
Matrix inline operator+(const double a, const Matrix& A ) {
  LOGInFct ( std::cout );
  size_t row =   A.getRow();
  size_t col =   A.getCol();
  double *tab = new double [ row * col ];
  const double *mat_ = A.getMatrix();

  for (size_t i=0; i < row*col; i++) {
    tab[i] = a + mat_[i];
  }
  Matrix tmp(row, col, tab, false);

  LOGOutFct( std::cout );
  return ( tmp );
}  

Matrix inline operator-(const double a, const Matrix& A ) {
  LOGInFct ( std::cout );
  size_t row =   A.getRow();
  size_t col =   A.getCol();
  double *tab = new double [ row * col ];
  const double *mat_ = A.getMatrix();

  for (size_t i=0; i < row*col; i++) {
    tab[i] = a - mat_[i];
  }
  Matrix tmp(row, col, tab, false);

  LOGOutFct( std::cout );
  return( tmp );
}

Matrix inline operator*(const double a, const Matrix& A ) {
  LOGInFct ( std::cout );
  size_t row =   A.getRow();
  size_t col =   A.getCol();
  double *tab = new double [ row * col ];
  const double *mat_ = A.getMatrix();

  for (size_t i=0; i < row*col; i++) {
    tab[i] = a * mat_[i];
  }
  Matrix tmp(row, col, tab, false);

  LOGOutFct( std::cout );
  return( tmp );
}

Matrix inline max( const double x, Matrix &A) {

  LOGInFct( std::cout );
  size_t row =   A.getRow();
  size_t col =   A.getCol();
  double *tab = new double [ row * col ];
  const double *mat_ = A.getMatrix();

  for (size_t i=0; i < row*col; i++) {
    tab[i] = MAX( mat_[i], x );
  }
  Matrix tmp(row, col, tab, false);

  LOGOutFct( std::cout );

  return ( tmp );
}

Matrix inline max( const Matrix &A, const double x) {

  LOGInFct( std::cout );

  size_t row =   A.getRow();
  size_t col =   A.getCol();
  double *tab = new double [ row * col ];
  const double *mat_ = A.getMatrix();

  for (size_t i=0; i < row*col; i++) {
    tab[i] = MAX( mat_[i], x );
  }
  Matrix tmp(row, col, tab, false);

  LOGOutFct( std::cout );

  return ( tmp );
}

Matrix inline max( const Matrix &A, Matrix &B ) {

  LOGInFct( std::cout );

  size_t row =   A.getRow();
  size_t col =   A.getCol();
  double *tab = new double [ row * col ];
  const double *matA_ = A.getMatrix();
  const double *matB_ = A.getMatrix();

  if( EQUAL_DIM_MAT( A, B ) ) {
    for (size_t i=0; i < row*col; i++) {
      tab[i] = MAX( matA_[i], matB_[i]);
    }
  } else {
    LOGMSG( 1, std::cout, not_same_dims, "" );
  }

  Matrix tmp(row, col, tab, false );
  LOGOutFct( std::cout );
  return ( tmp );
}


Matrix inline min(  const double x, const Matrix &A) {

  LOGInFct( std::cout );
  size_t row =   A.getRow();
  size_t col =   A.getCol();
  double *tab = new double [ row * col ];
  const double *mat_ = A.getMatrix();

  for (size_t i=0; i < row*col; i++) {
    tab[i] = MIN( mat_[i], x );
  }
  Matrix tmp(row, col, tab, false);

  LOGOutFct( std::cout );

  return ( tmp );
}

Matrix inline min( const Matrix &A, const double x) {

  LOGInFct( std::cout );

  size_t row =   A.getRow();
  size_t col =   A.getCol();
  double *tab = new double [ row * col ];
  const double *mat_ = A.getMatrix();

  for (size_t i=0; i < row*col; i++) {
    tab[i] = MIN( mat_[i], x );
  }
  Matrix tmp(row, col, tab, false);

  LOGOutFct( std::cout );

  return ( tmp );
}

Matrix inline min( const Matrix &A, const Matrix &B ) {

  LOGInFct( std::cout );

  size_t row =   A.getRow();
  size_t col =   A.getCol();
  double *tab = new double [ row * col ];
  const double *matA_ = A.getMatrix();
  const double *matB_ = A.getMatrix();

  if( EQUAL_DIM_MAT( A, B ) ) {
    for (size_t i=0; i < row*col; i++) {
      tab[i] = MIN( matA_[i], matB_[i]);
    }
  } else {
    LOGMSG( 1, std::cout, not_same_dims, "" );
  }

  Matrix tmp(row, col, tab, false);

  LOGOutFct( std::cout );

  return ( tmp );
}

double inline sum( const Matrix &A ) {

  size_t row =   A.getRow();
  size_t col =   A.getCol();
  double sum=0.0;

  const double *mat_ = A.getMatrix();

  for (size_t i=0; i < row*col; i++) {
    sum +=  mat_[i];
  }
  return( sum );
}

Matrix inline log( const Matrix &A ) {

  size_t row =   A.getRow();
  size_t col =   A.getCol();
  double *tab = new double [ row * col ];
  const double *mat_ = A.getMatrix();
  for (size_t i=0; i < row*col; i++) {
    tab[i] = log( mat_[i] );
  }
  Matrix tmp(row, col, tab, false);

  return( tmp );
}

Matrix inline abs( const Matrix &A ) {

  size_t row =   A.getRow();
  size_t col =   A.getCol();
  double *tab = new double [ row * col ];
  const double *mat_ = A.getMatrix();
  for (size_t i=0; i < row*col; i++) {
    tab[i] = fabs( mat_[i] );
  }
  Matrix tmp(row, col, tab, false);

  return( tmp );
}


#endif
