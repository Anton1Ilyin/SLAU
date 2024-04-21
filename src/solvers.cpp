#include"matrix.cpp"

std::vector<double> YacMult(CSR_matrix mat,std::vector<double> vec)
{
    std::vector<double> res(vec.size());
    for (int i=0; i<vec.size(); i++)
        {
            for(int j=mat.atRow(i);j<mat.atRow(i+1);j++)
            {
                res[i]+=mat.atVal(j)*vec[mat.atCol(j)];
            }
        }
    return res;
}
CSR_matrix invdiag(CSR_matrix A)
{
    int n=A.RowSize()-1;
    std::vector<double> diag(n);
    std::vector<unsigned> col(n), row(n + 1);
    for(int i=0; i<n; i++)
    {
        for(int j=A.atRow(i); j<A.atRow(i+1); j++)
        {
            if(i==A.atCol(j))
            {
                diag[i]=1/A.atVal(j);
                col[i]=i;
                row[i+1]=i+1;
            }
        }
    } 
}

std::vector<double> GZ(CSR_matrix& A, std::vector<double>& b, std::vector<double> x0, double breakpoint)
{
    std::vector<double> x=x0;
    CSR_matrix d=invdiag(A);
    double r=norm2(A*x-b);
    double z;
    while(breakpoint<r)
    {
        for(int i=0; i<x.size(); i++)
        {
            z=0.0;
            for(int j=A.atRow(i); j<A.atRow(i+1); j++)
            {
                if(i!=A.atCol(j))
                {
                    z+=A.atVal(j)*x[A.atCol(j)];
                }
            }
            x[i]=d.atVal(i)*(b[i]-z);
        }
        r=norm2(A*x-b);
    }
    return x;
};

std::vector<double> Yacoby(CSR_matrix& A, std::vector<double>& b, std::vector<double> x0, double breakpoint)
{
    std::vector<double> x=x0;
    double r=norm2(A*x-b);
    while(breakpoint<r)
    {
        x=(b-YacMult(A,x));
        r=norm2(A*x-b);
    }
    return x;
};

std::vector<double> Ineration(CSR_matrix& A, std::vector<double>& b, double tau, std::vector<double> x0, double breakpoint)
{
    std::vector<double> x=x0;
    std::vector<double> r=b-A*x;
    while(breakpoint<norm2(r))
    {
        x=x+r*tau;
        r=b-A*x;
    }
    return x;
};


std::vector<double> Chebyshev(CSR_matrix& A, std::vector<double>& b, double lambda_min, double lambda_max, std::vector<double> x0, double breakpoint)
{
    std::vector<double> x=x0;
    std::vector<double> r=b-A*x;
    int Tn=128;
    int deg=7;
    double c=cos(M_PI/Tn);
    double s=sin(M_PI/Tn);
    std::vector<double> root(Tn);
    std::vector<int> index(Tn);
    index[0]=0;
    int k=Tn/4;
    index[Tn/2]=1;
    for(int i=2; i<deg; i++)
    {
        for(int j=0; j<Tn; j+=k)
        {
            index[j+k/2]=Tn/(k/2)-1-index[j];
        }
        k/=2;
    }
    root[0]=(lambda_max+lambda_min)/2+(lambda_max-lambda_min)*cos(M_PI/(2*Tn))/2;
    for(int i=1; i<Tn; i++)
    {
        root[i]=(lambda_max+lambda_min)/2+(lambda_max-lambda_min)/2*(root[i-1]*c-s*sqrt(1-root[i-1]*root[i]));
    }
    while(breakpoint<norm2(r))
    {
        for(int i=0; i<Tn; i++)
        {
            x=x+r/root[index[i]];
            r=b-A*x;
        }
    }
    return x;
};

std::vector<double> SGZ(CSR_matrix& A, std::vector<double>& b, std::vector<double> x0, double breakpoint)
{
    std::vector<double> x=x0;
    CSR_matrix d=invdiag(A);
    double r=norm2(A*x-b);
    double z;
    while(breakpoint<r)
    {
        for(int i=0; i<x.size(); i++)
        {
            z=0;
            for(int j=A.atRow(i); j<A.atRow(i+1); j++)
            {
                if(i!=A.atCol(j))
                {
                    z+=A.atVal(j)*x[A.atCol(j)];
                }
            }
            x[i]=d.atVal(i)*(b[i]-z);
        }
        for(int i=x.size()-1; i>0; i--)
        {
            z=0;
            for(int j=A.atRow(i); j<A.atRow(i+1); j++)
            {
                if(i!=A.atCol(j))
                {
                    z+=A.atVal(j)*x[A.atCol(j)];
                }
            }
            x[i]=d.atVal(i)*(b[i]-z);
        }
        r=norm2(A*x-b);
    }
    return x;
};

std::vector<double> Gradient(CSR_matrix& A, std::vector<double>& b, std::vector<double> x0, double breakpoint)
{
    std::vector<double> x=x0;
    std::vector<double> r=b-A*x;
    double tau=r*r/(r*(A*r));
    while(breakpoint<norm2(r))
    {
        x=x+r*tau;
        r=b-A*x;
        tau=r*r/(r*(A*r));
    }
    return x;
};

std::vector<double> ChebSGZ(CSR_matrix& A, std::vector<double>& b, int n, double rho, std::vector<double> x0, double breakpoint)
{
    std::vector<double> x=x0;
    CSR_matrix d=invdiag(A);
    double r=norm2(A*x-b);
    double z;
    std::vector<double> mu(n);
    mu[0]=1;
    mu[1]=1/rho;
    for(int i=2; i<n; i++)
    {
        mu[i]=2/rho*(mu[i-1]-mu[i-2]);
    }
    std::vector<double> y0=x0;
    std::vector<double> y1;
    for(int i=0; i<x.size(); i++)
    {
        z=0;
        for(int j=A.atRow(i); j<A.atRow(i+1); j++)
        {
            if(i!=A.atCol(j))
            {
                z+=A.atVal(j)*x[A.atCol(j)];
            }
        }
        x[i]=d.atVal(i)*(b[i]-z);
    }
    for(int i=x.size()-1; i>0; i--)
    {
        z=0;
        for(int j=A.atRow(i); j<A.atRow(i+1); j++)
        {
            if(i!=A.atCol(j))
            {
                z+=A.atVal(j)*x[A.atCol(j)];
            }
        }
        x[i]=d.atVal(i)*(b[i]-z);
    }
    y1=x;
    r=norm2(A*x-b);
    while(breakpoint<r)
    {
        for(int k=0; k<n-1; k++)
        {
            for(int i=0; i<x.size(); i++)
            {
                z=0;
                for(int j=A.atRow(i); j<A.atRow(i+1); j++)
                {
                    if(i!=A.atCol(j))
                    {
                        z+=A.atVal(j)*x[A.atCol(j)];
                    }
                }
                x[i]=d.atVal(i)*(b[i]-z);
            }
            for(int i=x.size()-1; i>0; i--)
            {
                z=0;
                for(int j=A.atRow(i); j<A.atRow(i+1); j++)
                {
                    if(i!=A.atCol(j))
                    {
                        z+=A.atVal(j)*x[A.atCol(j)];
                    }
                }
                x[i]=d.atVal(i)*(b[i]-z);
            }
            x=x*(2*mu[k]/rho/mu[k+1])-y0*(mu[k-1]/mu[k+1]/rho);
            y0=y1;
            y1=x;
            r=norm2(A*x-b);
        } 
    }
    return x;
};

std::vector<double> ConjGradient(CSR_matrix& A, std::vector<double>& b, std::vector<double> x0, double breakpoint)
{
    std::vector<double> x=x0;
    std::vector<double> r=b-A*x;
    std::vector<double> d=r;
    double tau;
    while(breakpoint<norm2(r))
    {
        tau=r*r/(d*(A*d));
        x=x-d*tau;
        r=b-A*x;
        if(norm2(d)==0)
        {
            break;
        }
        else
        {
            tau=d*d/(r*r);
            d=r+d*tau;
        }
    }
    return x;
};

class Arnoldi
{
    private:
    std::vector<std::vector<double>> vecs;
    CSR_matrix H;
    public:
    CSR_matrix get_H()
    {
        return H;
    }
    std::vector<std::vector<double>> get_vecs()
    {
        return vecs;
    }
    std::vector<double> get_vec(int i)
    {
        return vecs[i];
    }
    Arnoldi(CSR_matrix A, std::vector<double> b, std::vector<double> x0, int i)
    {
        std::vector<double> v=(A*x0-b)/norm2(A*x0-b);
        std::vector<double> t;
        std::vector<int> col, row;
        std::vector<double> val;
        vecs.push_back(v);
        row.push_back(0);
        int c=0;
        for(int j=0; j<i; j++)
        {
            t=A*v;
            double h;
            for(int k=0; k<=j; k++)
            {
                h=v*t;
                val.push_back(h);
                col.push_back(k);
                c++;
                t=t-v*h;
            }
            val.push_back(norm2(t));
            col.push_back(j+1);
            c++;
            row.push_back(c);
            v=t/h;
            vecs.push_back(v);
        }
        H=CSR_matrix(val,col,row).transpose();
    }
};