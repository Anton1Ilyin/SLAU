#include<vector>
#include<cmath>

class vector_ND
{
    private:
    std::vector<double> data;
    public:
    vector_ND(std::vector<double> v)
    {
        data=v;
    }
    std::size_t size()
    {
        return data.size();
    }
    std::vector<double> get_data()
    {
        return data;
    }
    const vector_ND& operator+(const vector_ND& vec)
    {
        std::vector<double> res(data.size());
        for(int i=0;i<data.size();i++)
        {
            res[i]=data[i]+vec.data[i];
        }
        return vector_ND(res);
    }
    const vector_ND& operator-(const vector_ND& vec)
    {
        std::vector<double> res(data.size());
        for(int i=0;i<data.size();i++)
        {
            res[i]=data[i]-vec.data[i];
        }
        return vector_ND(res);
    }
    double operator*(const vector_ND& vec) const
    {
        double res;
        for(int i=0; i<vec.data.size();i++)
        {
            res+=data[i]+vec.data[i];
        }
        return res;
    }
    const vector_ND& operator*(double a)
    {
        std::vector<double> res(data.size());
        for(int i=0;i<data.size();i++)
        {
           res[i]=a*data[i];
        }
        return vector_ND(res);
    }
};

class matrix
{
    std::vector<double> data;
    int col;
    
    public:
    const double operator()(std::size_t i, std::size_t j) const {
        return data[i*col+j];
    }
    double get_col()
    {
        return col;
    }
    unsigned int data_size()
    {
        return data.size();
    }
    matrix(unsigned int i, unsigned int j)
    {
        std::vector<double> data(i*j);
        col=j;
    }
    void fill(std::vector<double> d)
    {
        data=d;
    }
    const matrix& operator+(const matrix& mat)
    {
        matrix res(data.size()/col, col);
        std::vector<double> d(data.size());
        for(int i=0;i<data.size();i++)
        {
            d[i]=data[i]+mat.data[i];
        }
        res.fill(d);
        return res;
    }
    const matrix& operator-(const matrix& mat)
    {
        matrix res(data.size()/col, col);
        std::vector<double> d(data.size());
        for(int i=0;i<data.size();i++)
        {
            d[i]=data[i]+mat.data[i];
        }
        res.fill(d);
        return res;
    }
    std::vector<double> operator*(const std::vector<double>& vec) const
    {
        std::vector<double> res(vec.size());
        for(int i=0; i<vec.size();i++)
        {
            double a;
            for(int j=0;j<vec.size();j++)
            {
                a+=vec[j]*data[i*col+j];
            }
            res[i]=a;
        }
        return res;
    }
    // vector_ND operator*(const vector_ND& vec) const
    // {
    //     std::vector<double> res(vec.size());
    //     for(unsigned int i=0; i<vec.size();i++)
    //     {
    //         double a;
    //         for(int j=0;j<vec.size();j++)
    //         {
    //             a+=vec.data[j]*data[i*col+j];
    //         }
    //         res[i]=a;
    //     }
    //     return vector_ND(res);
    // }
    const matrix& operator*(double a)
    {
        matrix res(data.size()/col, col);
        std::vector<double> d(data.size());
        for(int i=0;i<data.size();i++)
        {
            d[i]=data[i]*a;
        }
        res.fill(d);
        return res;
    }
    const matrix& operator*(const matrix& mat)
    {
        matrix res(col, mat.data.size()/mat.col);
        std::vector<double> d(data.size());
        for(int i=0;i<col;i++)
        {
            for(int j=0; j<mat.data.size()/mat.col; j++)
            {
                double a=0;
                for(int k=0; k<col; k++)
                {
                    a+=data[i*col+k]*mat(k,j);
                }
                d[i*col+j]=a;
            }
        }
        res.fill(d);
        return res;
    }
    matrix transpose()
    {
        matrix res(data.size()/col, col);
        std::vector<double> d(data.size());
        for(int i=0;i<data.size()/col;i++)
        {
            for(int j=0;j<col;j++)
            {
                d[j*col+i]=data[i*col+j];
            }
        }
        res.fill(d);
        return res;
    }
    void ones(std::size_t n)
    {
        data.resize(n*n);
        col=n;
        for(std::size_t i=0; i<n; i++)
        {
            for(std::size_t j=0; j<n; j++)
            {
                if (i!=j)
                {
                    data[i*n+j]=0;
                }
                else
                {
                    data[i*n+i]=1;
                }
            }
        }
    }
};

class CSR_matrix
{
    private:
    std::vector<double> val;
    std::vector<int> column;
    std::vector<int> row;
    public:
    int RowSize()
    {
        return row.size();
    }
    double atVal(int i)
    {
        return val[i];
    }
    int atCol(int i)
    {
        return column[i];
    }
    int atRow(int i)
    {
        return row[i];
    }
    CSR_matrix(matrix data)
    {
        this->row.push_back(0);
        int k=0;
        for(int i=0;i<data.data_size()/data.get_col();i++)
        {
            
            for(int j=0; j<data.get_col();j++)
            {
                double d=data(i,j);
                if (d!=0)
                {
                    this->val.push_back(d);
                    this->column.push_back(j);
                    k++;
                }
            }
            this->row.push_back(k);
        }
    }
    CSR_matrix(std::vector<double> v, std::vector<int> c, std::vector<int> r)
    {
        val=v;
        column=c;
        row=r;
    }
    double operator()(int i, int j) const
    {
        for(int k=row[i]; k<row[i+1]; k++)
        {
            if(column[k]==j)
            {
                return val[k];
            }
        }
        return 0;
    }
    std::vector<double> operator*(const std::vector<double>& vec) const
    {
        std::vector<double> res(vec.size());
        for (int i=0; i<vec.size(); i++)
        {
            for(int j=row[i];j<row[i+1];j++)
            {
                res[i]+=val[j]*vec[column[j]];
            }
        }
        return res;
    }
    // vector_ND operator*(const vector_ND& vec) const
    // {
    //     std::vector<double> res(vec.data.size());
    //     for(int i=0; i<vec.data.size();i++)
    //     {
    //         double a;
    //         for(int j=0;j<vec.data.size())
    //         {
    //             a+=vec.data[j]*(*this(i,j));
    //         }
    //         res[i]=a;
    //     }
    //     return vector_ND(res);
    // }
    const CSR_matrix& operator*(double a)
    {
        std::vector<double> res;
        for(int i=0; i<val.size();i++)
        {
            res.push_back(val[i]*a);
        }
        return CSR_matrix(res,column,row);
    }
    const CSR_matrix& operator*(const CSR_matrix& mat)
    {
        std::vector<double> v;
        std::vector<int> c;
        std::vector<int> r;
        r.push_back(0);
        int l;
        double res;
        for (int i=0; i<row.size()-1; i++)
        {
            l=0;
            for(int k=mat.row[i];k<mat.row[i+1];k++)
            {
                res=0.0;
                for(int j=row[i];j<row[i+1];j++)
                {
                    if(column[k]==j)
                    {
                        res+=val[k]*mat(i,j);
                    }
                }
                if(res!=0.0)
                {
                    v.push_back(res);
                    c.push_back(k);
                    l++;
                }
            }
            r.push_back(l);
        }
        return CSR_matrix(v,c,r);
    }
    CSR_matrix transpose()
    {
        int N=0;
        for(int i=0; i<column.size(); i++)
        if(column[i]>N)
        N=column[i];
        N++;
        std::vector<double> vT; 
        std::vector<int> cT, rT;
        for(int i=0; i<val.size(); i++)
        {
            if(column[i]+2<N+1)
            {
                rT[column[i]+2]++;
            }
        }
        for(int i=1;i<N;i++)
        {
            rT[i+1]+=rT[i];
        }
        for(int i=0; i<row.size()-1;i++)
        {
            for(int j=row[i]; j<row[i+1]; j++)
            {
                int k=rT[column[i]+1];
                vT[k]=val[i];
                cT[k]=i;
                rT[column[i]+1]++;
            }
        }
        return CSR_matrix(vT, cT, rT);
    }
    CSR_matrix()
    {
        row={};
        val={};
        column={};
    } 
    CSR_matrix& operator=(const CSR_matrix& right)
    {
        val=right.val;
        column=right.column;
        row=right.row;
    }
    CSR_matrix(int j, int i, double alpha, double beta)
    {
        double s=beta/sqrt(alpha*alpha+beta*beta), c=alpha/sqrt(alpha*alpha+beta*beta);
        row.push_back(0);
        for(int k=0; k<i; k++)
        {
            if(k==j)
            {
                val.push_back(c);
                column.push_back(k);
                val.push_back(s);
                column.push_back(k+1);
                row.push_back(row[k]+2);
            }
            else if(k==j+1)
            {
                val.push_back(-s);
                column.push_back(k-1);
                val.push_back(c);
                column.push_back(k);
                row.push_back(row[k]+2);
            }
            else
            {
                val.push_back(1);
                column.push_back(k);
                row.push_back(row[k]+1);
            }
        }
    }
};




std::vector<double> operator+(std::vector<double> left, std::vector<double> right)
{
    std::vector<double> res(left.size());
    for(int i=0; i<left.size(); i++)
    {
        res[i]=left[i]+right[i];
    }
    return res;
}
std::vector<double> operator-(std::vector<double> left, std::vector<double> right)
{
    std::vector<double> res(left.size());
    for(int i=0; i<left.size(); i++)
    {
        res[i]=left[i]-right[i];
    }
    return res;
}
double operator*(const std::vector<double>& left, const std::vector<double>& right) {
    double res;
    for (unsigned int i=0; i<left.size(); i++)
    {
        res+=left[i]+right[i];
    }
    return res;
};
std::vector<double> operator*(const std::vector<double>& vec, double a) {
    unsigned n = vec.size();
    std::vector<double> res(n);

    for (unsigned i = 0; i < vec.size(); i++) {
        res[i] = vec[i] * a;
    }
    return res;
}
std::vector<double> operator/(const std::vector<double>& vec, double a) {
    unsigned n = vec.size();
    std::vector<double> res(n);

    for (unsigned i = 0; i < vec.size(); i++) {
        res[i] = vec[i] / a;
    }
    return res;
}
double norm2(const std::vector<double>& vec) {
    return sqrt(vec * vec);
}