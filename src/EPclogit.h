class EPlogit
{
	public:
		EPlogit(void)
		{}
		void operator()(mat X,mat Yt, double damp,double lambda)
		{
			int p=X.n_cols;
			int n=X.n_rows;
			mat mu(p,1);
			mat Q(p,p);
			mat I(p,p);
			I.eye();
			for(int i=0;i<n;i++)
			{
				if(Yt(i,0)==0)
				{
					X(i,span::all)=-X(i,span::all);
				}
			}
			Q.eye();
			Q=(double)1/lambda*Q;
			mat K(n,1);
			K.fill(0);
			mat H(n,1);
			H.fill(0);
			mat m(p,1);
			mat Qm(p,1);
			m.fill(0);
			Qm.fill(0);
			mat V=lambda*I;
			int count=0;
			mat Qmt=Qm;
			double crit=1;
			while(crit>0.01 & count< 10)
			{
				_Z=0;
				for(int i=0;i<n;i++)
				{
					mat x=X(i,span::all);
					double k=K(i,0);	
					double h=H(i,0);
					mat Qmumi=Qm-x.t()*h;
					mat Qmi=Q-k*x.t()*x;
					mat Vmi=V+k*(V*x.t()*x*V)/as_scalar((1-k*x*V*x.t()));
					mat mmi=Vmi*Qmumi;
					double  mi=dot(mmi,x);
					double  vi=dot(x*Vmi,x);
					double z=1;
					if((vi>-1)){
						double ratio=mi/sqrt((PI/8)*vi+1);
						double t;
						z=logit(ratio);
						t=dlogit(ratio);
						if(z<0.0000001) z=0.0000001;
						if(z>0.9999999) z=0.9999999;
						double dzm=t/sqrt((PI/8)*vi+1);
						double dzv=-0.5*(PI/8)*t*ratio/((PI/8)*vi+1);
						double mn=mi+vi*dzm/z;
						double vn=vi-vi*(dzm*dzm/(z*z)-2*dzv/z)*vi;
						K(i,0)=damp*(1.0/vn-1.0/vi)+(1-damp)*k;
						H(i,0)=damp*(mn/vn-mi/vi)+(1-damp)*h;
						
						k=K(i,0);	
						h=H(i,0);
						double Kmi=(double)1/vi;
						double hmi=vi/mi;
						V=Vmi-Vmi*x.t()*x*Vmi/(1.0/K(i,0)+as_scalar(x*Vmi*x.t()));
						mat vv=x.t()*x*K(i,0);
						mat mm=H(i,0)*x.t();
						Q=Qmi+x.t()*x*K(i,0);
						Qm=Qmumi+H(i,0)*x.t();
						m=V*Qm;
						double r=mi/vi;
						double q=1.0/vi;
						_Z+=log(z)+0.5*mi*mi/(vi)-0.5*(h+r)*(h+r)/(k+q)+0.5*log(vi*(k+q));
					}
				}
				count++;
				crit=norm(Qm-Qmt)/max(norm(Qmt),max(norm(Qm),0.001));
				Qmt=Qm;
				double val1,s1;
				log_det(val1,s1,V);
				_Z+=0.5*val1-p*0.5*log(lambda)+as_scalar(0.5*m.t()*Q*m);
			

			
			}
			_m=m;
			_V=V;

		}
		double logit(double x)
		{
			return 1.0/(1+exp(-x));
		}
		double dlogit(double x)
		{
			return exp(-x)/((1+exp(-x))*(1+exp(-x)));
		}
		mat Get_m(void){return _m;}
		mat Get_V(void){return _V;}
		double Get_Z(void){return _Z;}

	private:
		mat _m;
		boost::math::normal_distribution<> *_d;
		mat _V;
		double _Z;
};
