#ifndef GLM_H
#define GLM_H

/**
 * @brief The GLM class
 * Generalized linear model
 */
class GLM
{
    public:
        /**
         * @brief linkfunc
         *        link function g(u)=Xb
         * @param u
         * @return
         */
        static double linkfunc(double u);
        /**
         * @brief dlinkfunc
         *        first order derivative of link function
         * @param u
         * @return
         */
        static double dlinkfunc(double u);
        /**
         * @brief ddlinkfunc
         *        second order derivative of link function
         * @param u
         * @return
         */
        static double ddlinkfunc(double u);
        /**
         * @brief invlinkfunc
         *        inverse of link function
         * @param Xb
         * @return
         */
        static double invlinkfunc(double Xb);

    public:
        /**
         * @brief irls
         * @param X
         * @param y
         * @param M
         * @param N
         * @param b
         * @param maxIter
         * @param epsilon
         */
        static void irls(const double* X, const double* y, int N,
                         double *b, int maxIter, double epsilon);
        /**
         * @brief irls
         * @param X
         * @param y
         * @param N
         * @param b0
         * @param b
         * @param maxIter
         * @param epsilon
         */
        static void irls(const double* X, const double* y, int N,
                         const double* b0, double* b,
                         int maxIter, double epsilon);
        /**
         * @brief sse
         * @param X
         * @param y
         * @param N
         * @param b
         * @return
         */
        static double sse(const double* X, const double* y, int N, const double* b);

    public:
        GLM();
        virtual ~GLM();
};

#endif // GLM_H
