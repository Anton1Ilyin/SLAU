#include<vector>

class vector_ND
{
    private:
    std::vector<double> data;
    public:
    vector_ND(std::vector<double> v)
    {
        data=v;
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
        for(int i=0; i<vec.size();i++)
        {
            res+=data[i]+vec.data[i]
        }
        return res;
    }
    const vector_ND& operator*(double a)
    {
        std::vector res(data.size());
        for(int i=0;i<data.size();i++)
        {
           res[i]=a*data[i];
        }
        return vector_ND(res);
    }
}

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
    const matrix& operator+(const matrix& mat)
    {
        matrix res(data.size()/col, col);
        for(int i=0;i<data.size()/col;i++)
        {
            for(int j=0;j<col;j++)
            {
                res(i,j)=data[i*col+j]+mat(i,j);
            }
        }
        return res;
    }
    const matrix& operator-(const matrix& mat)
    {
        matrix res(data.size()/col, col);
        for(int i=0;i<data.size()/col;i++)
        {
            for(int j=0;j<col;j++)
            {
                res(i,j)=data[i*col+j]-mat(i,j);
            }
        }
        return res;
    }
    std::vector<double> operator*(const std::vector<double>& vec) const
    {
        std::vector<double> res(vec.size());
        for(int i=0; i<vec.size();i++)
        {
            double a;
            for(int j=0;j<vec.size())
            {
                a+=vec[j]*data[i*col+j];
            }
            res[i]=a;
        }
        return res;
    }
    vector_ND operator*(const vector_ND& vec) const
    {
        std::vector<double> res(vec.data.size());
        for(int i=0; i<vec.data.size();i++)
        {
            double a;
            for(int j=0;j<vec.data.size())
            {
                a+=vec.data[j]*data[i*col+j];
            }
            res[i]=a;
        }
        return vector_ND(res);
    }
    const matrix& operator*(double a)
    {
        matrix res(data.size()/col, col);
        for(int i=0;i<data.size()/col;i++)
        {
            for (int j=0; j<col; j++)
            {
                res(i,j)=data[i*col+j]*a;
            }
        }
        return res;
    }
    const matrix& operator*(const matrix& mat)
    {
        matrix res(col, mat.data.size()/mat.col);
        for(int i=0;i<col;i++)
        {
            for(int j=0; j<mat.data.size()/mat.col; j++)
            {
                double a=0;
                for(int k=0; k<col; k++)
                {
                    a+=data[i*col+k]*mat(k,j);
                }
                res(i,j)=a;
            }
        }
        return res;
    }
}

class CSR_matrix
{
    private:
    std::vector<double> val;
    std::vector<int> column;
    std::vector<int> row;
    public:
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
    const double operator()(std::size_t i, std::size_t j) const
    {
        if(row[i]==row[i+1])
        {
            return 0.0;
        }
        else if (column[row[i]+j]==j)
        {
            return val[row[i]+j]
        }
        else
        {
            return 0.0;
        }
    }
    std::vector<double> operator*(const std::vector<double>& vec) const
    {
        std::vector<double> res(vec.size());
        for(int i=0; i<vec.size();i++)
        {
            double a;
            for(int j=0;j<vec.size())
            {
                a+=vec[j]*(*this(i,j));
            }
            res[i]=a;
        }
        return res;
    }
    vector_ND operator*(const vector_ND& vec) const
    {
        std::vector<double> res(vec.data.size());
        for(int i=0; i<vec.data.size();i++)
        {
            double a;
            for(int j=0;j<vec.data.size())
            {
                a+=vec.data[j]*(*this(i,j));
            }
            res[i]=a;
        }
        return vector_ND(res);
    }
    const matrix& operator*(double a)
    {
        std::vector<double> res;
        for(int i=0; i<val.size();i++)
        {
            res.push_back(val[i]*a);
        }
        return CSR_matrix(res,column,row);
    }
}
