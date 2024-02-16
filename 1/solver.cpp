#include <fstream>
#include <vector>


std::vector<double> solver(std::string name)
{
    std::ifstream file;
	file.open(name);
    int N;
    file>>N;
    std::vector<double> a(N);
    std::vector<double> b(N+1);
    std::vector<double> c(N);
    std::vector<double> d(N+1);
    for(int i=0;i<N;i++)
    {
        file>>a[i];
    }
    for(int i=0;i<N+1;i++)
    {
        file>>b[i];
    }
    for(int i=0;i<N;i++)
    {
        file>>c[i];
    }
    for(int i=0;i<N+1;i++)
    {
        file>>d[i];
    }
    std::vector<double> p(N);
    std::vector<double> q(N);
    std::vector<double> x(N+1);
    p[0]=-c[0]/b[0];
    q[0]=d[0]/b[0];
    for(int i=1; i<N; i++)
    {
        p[i]=-c[i]/(a[i-1]*p[i-1]+b[i]);
        q[i]=(d[i]-a[i-1]*q[i-1])/(a[i]*p[i-1]+b[i]);
    }
    x[N]=(d[N]-a[N-1]*q[N-1])/(a[N-1]*p[N-1]+b[N]);
    for(int i=N-1;i>=0;i--)
    {
        x[i]=p[i]*x[i+1]+q[i];
    }
    return x;
}
