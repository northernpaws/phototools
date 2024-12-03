

#ifndef PHOTOTOOLS_MATRIX_H
#define PHOTOTOOLS_MATRIX_H

#include <type_traits>
#include <stdexcept>
#include <array>
#include <initializer_list>
#include "../color.h"

namespace core {

    template<class T, size_t _X, size_t _Y>
    struct Matrix {
        static_assert(std::is_arithmetic<T>::value, "T must be numeric");
        static_assert(std::is_floating_point<T>::value, "T must be a floating point type");

        static_assert(_X > 1, "_X must be greater then 1");
        static_assert(_Y > 1, "_Y must be greater then 1");

        typedef std::array<T, _X> Row;

        // reasons for using a C++ array over C-style: https://stackoverflow.com/a/77106115
        std::array<Row, _Y> _rows;

        Matrix() = default;

        /**
         * Constructor for creating matrices from C-style arrays.
         * @param mat
         */
        Matrix(T mat[_Y][_X]) {
            for (int y = 0; y < _Y; y++) {
                for (int x = 0; x < _X; x++) {
                    _rows[y][x] = mat[y][x];
                }
            }
        }

        /**
         * Construct for creating matrixes directly from C++ arrays.
         *
         * This also allows for the use of initialization list syntax, for example:
         *  Matrix3x3d({ {1, 1, 1}, {1, 1, 1}, {1, 1, 1} })
         *
         * @param mat
         */
        Matrix(std::array<Row, _Y> mat) {
            // TODO: does this need a copy operator?
            _rows = mat;
        }

        /**
         * Copy constructor.
         * @param mat
         */
        Matrix(const Matrix &mat) {
            // TODO: ensure this is doing a copy
            _rows = mat._rows;
        }

        /**
         * Constructor for using initializer lists, i.e.:
         * Matrix3x3d({ {1, 1, 1}, {1, 1, 1}, {1, 1, 1} })
         *
         * @param lst
         */
        // TODO: may be able to remove in favor of the above std::array constructor which should allow the same syntax.
        Matrix(std::initializer_list<std::initializer_list<T>> lst) {
            assert(lst.size() == _Y);
//        static_assert(lst.size() > _Y, "Initializer list is taller then _Y");
//        static_assert(lst[0].size() > _X, "Initializer list is wider then _X");

            int y = 0;
            for (const std::initializer_list<T> &row: lst) {
                assert(row.size() == _X);

                std::copy(row.begin(), row.end(), _rows[y].begin());

                y += 1;
            }
        }

        Row &operator[](unsigned int index) {
            return _rows[index];
        }

        Row operator[](unsigned int index) const {
            return _rows[index];
        }

        /**
         * Assign from C-style array.
         *
         * @param mat
         * @return
         */
        Matrix &operator=(T mat[_Y][_X]) {
            for (int y = 0; y < _Y; y++) {
                for (int x = 0; x < _X; x++) {
                    _rows[y][x] = mat[y][x];
                }
            }
            return *this;
        }
    };

// define templates for common matrix dimensions

    template<typename T>
    struct Matrix3x3 : public Matrix<T, 3, 3> {
        using typename Matrix<T, 3, 3>::Row;

        Matrix3x3() = default;

        /**
         * Convenience constructor for initializer list syntax.
         */
        Matrix3x3(Row row1, Row row2, Row row3) {
            this->_rows[0] = row1;
            this->_rows[1] = row2;
            this->_rows[2] = row3;
        }

        inline Color transform(float r, float g, float b) const {
            Color transformed;

            transformed.r = this->_rows[0][0] * r + this->_rows[0][1] * g + this->_rows[0][2] * b;
            transformed.r = this->_rows[1][0] * r + this->_rows[1][1] * g + this->_rows[1][2] * b;
            transformed.r = this->_rows[2][0] * r + this->_rows[2][1] * g + this->_rows[2][2] * b;

            return transformed;
        }

        /**
         * Transforms a color using the matrix where each matrix
         * row is a color component and each column is a pixel.
         *
         * @param c
         * @return
         */
        inline Color transform(const Color &c) const {
            Color transformed;

            transformed.r = this->_rows[0][0] * c.r + this->_rows[0][1] * c.g + this->_rows[0][2] * c.b;
            transformed.r = this->_rows[1][0] * c.r + this->_rows[1][1] * c.g + this->_rows[1][2] * c.b;
            transformed.r = this->_rows[2][0] * c.r + this->_rows[2][1] * c.g + this->_rows[2][2] * c.b;

            return transformed;
        }

        // ref: https://stackoverflow.com/a/18504573
// TODO: original has the note: "The computations are correct, but result become unstable when determinant is close to
//  0. I found this by comparing it with with the numpy.linalg.inv method which provides much better result. Just
//  multiply the inverted and the original matrices and compare it with the identity matrix. "
        inline Matrix3x3<T> inverse() const {
            auto _rows = Matrix<T, 3, 3>::_rows;

            T determinant = _rows[0][0] * (_rows[1][1] * _rows[2][2] - _rows[2][1] * _rows[1][2]) -
                            _rows[0][1] * (_rows[1][0] * _rows[2][2] - _rows[1][2] * _rows[2][0]) +
                            _rows[0][2] * (_rows[1][0] * _rows[2][1] - _rows[1][1] * _rows[2][0]);

            if (determinant == 0) {
                // TODO: document in here why this is a problem, re: https://stackoverflow.com/a/18504573 (maybe https://stackoverflow.com/a/984286)
                throw std::runtime_error("inverse determinant is 0");
            }

            T inverse_determinant = 1 / determinant;

            Matrix3x3<T> inverted;

            inverted[0][0] = (_rows[1][1] * _rows[2][2] - _rows[2][1] * _rows[1][2]) * inverse_determinant;
            inverted[0][1] = (_rows[0][2] * _rows[2][1] - _rows[0][1] * _rows[2][2]) * inverse_determinant;
            inverted[0][2] = (_rows[0][1] * _rows[1][2] - _rows[0][2] * _rows[1][1]) * inverse_determinant;
            inverted[1][0] = (_rows[1][2] * _rows[2][0] - _rows[1][0] * _rows[2][2]) * inverse_determinant;
            inverted[1][1] = (_rows[0][0] * _rows[2][2] - _rows[0][2] * _rows[2][0]) * inverse_determinant;
            inverted[1][2] = (_rows[1][0] * _rows[0][2] - _rows[0][0] * _rows[1][2]) * inverse_determinant;
            inverted[2][0] = (_rows[1][0] * _rows[2][1] - _rows[2][0] * _rows[1][1]) * inverse_determinant;
            inverted[2][1] = (_rows[2][0] * _rows[0][1] - _rows[0][0] * _rows[2][1]) * inverse_determinant;
            inverted[2][2] = (_rows[0][0] * _rows[1][1] - _rows[1][0] * _rows[0][1]) * inverse_determinant;

            return inverted;
        }
    };

    template<typename T>
    struct Matrix4x3 : public Matrix<T, 4, 3> {
    };

// define common matrix types

    typedef Matrix3x3<float> Matrix3x3f;
    typedef Matrix3x3<double> Matrix3x3d;

    typedef Matrix4x3<float> Matrix4x3f;
    typedef Matrix4x3<double> Matrix4x3d;

// ref: https://stackoverflow.com/a/18504573
// TODO: original has the note: "The computations are correct, but result become unstable when determinant is close to
//  0. I found this by comparing it with with the numpy.linalg.inv method which provides much better result. Just
//  multiply the inverted and the original matrices and compare it with the identity matrix. "
    template<typename T>
    Matrix3x3<T> inverse(const Matrix3x3<T> &mat) {
        T determinant = mat[0][0] * (mat[1][1] * mat[2][2] - mat[2][1] * mat[1][2]) -
                        mat[0][1] * (mat[1][0] * mat[2][2] - mat[1][2] * mat[2][0]) +
                        mat[0][2] * (mat[1][0] * mat[2][1] - mat[1][1] * mat[2][0]);

        if (determinant == 0) {
            // TODO: document in here why this is a problem, re: https://stackoverflow.com/a/18504573 (maybe https://stackoverflow.com/a/984286)
            throw std::runtime_error("inverse determinant is 0");
        }

        T inverse_determinant = 1 / determinant;

        Matrix3x3<T> inverted;
        inverted[0][0] = (mat[1][1] * mat[2][2] - mat[2][1] * mat[1][2]) * inverse_determinant;
        inverted[0][1] = (mat[0][2] * mat[2][1] - mat[0][1] * mat[2][2]) * inverse_determinant;
        inverted[0][2] = (mat[0][1] * mat[1][2] - mat[0][2] * mat[1][1]) * inverse_determinant;
        inverted[1][0] = (mat[1][2] * mat[2][0] - mat[1][0] * mat[2][2]) * inverse_determinant;
        inverted[1][1] = (mat[0][0] * mat[2][2] - mat[0][2] * mat[2][0]) * inverse_determinant;
        inverted[1][2] = (mat[1][0] * mat[0][2] - mat[0][0] * mat[1][2]) * inverse_determinant;
        inverted[2][0] = (mat[1][0] * mat[2][1] - mat[2][0] * mat[1][1]) * inverse_determinant;
        inverted[2][1] = (mat[2][0] * mat[0][1] - mat[0][0] * mat[2][1]) * inverse_determinant;
        inverted[2][2] = (mat[0][0] * mat[1][1] - mat[1][0] * mat[0][1]) * inverse_determinant;

        return inverted;
    }
}

#endif //PHOTOTOOLS_MATRIX_H
