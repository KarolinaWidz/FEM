import static java.lang.Math.sqrt;

class UniversalElement {
	private IntegralPoint [] integralPoints;
	private int numberOfIntegralPoints;
	private double [][] dNdksi;
	private double [][] dNdeta;


	UniversalElement() {
		this.integralPoints = setIntegralPoints();
		this.numberOfIntegralPoints = 4;
		this.dNdksi= setDNdksi();
		this.dNdeta = setDNdeta();

	}

	private IntegralPoint[] setIntegralPoints() { //Add integral points to universal element
		IntegralPoint [] integralPoints = new IntegralPoint[4];
		integralPoints[0] = new IntegralPoint(-1/sqrt(3),-1/sqrt(3),1,1);
		integralPoints[1] = new IntegralPoint(1/sqrt(3),-1/sqrt(3),1,1);
		integralPoints[2] = new IntegralPoint(1/sqrt(3),1/sqrt(3),1,1);
		integralPoints[3] = new IntegralPoint(-1/sqrt(3),1/sqrt(3),1,1);
		return integralPoints;
	}

	private double [] setShapeFunctions(double ksi, double eta) { //set local shape functions
		double [] tmp = new double[4];
		tmp[0] = ((1-ksi)*(1-eta))/4;
		tmp[1] = ((1+ksi)*(1-eta))/4 ;
		tmp[2] = ((1+ksi)*(1+eta))/4 ;
		tmp[3] = ((1-ksi)*(1+eta))/4 ;
		return tmp;
	}

	void setShapeFunctionInIntegralPoints(){ // return array with shape function in integral points
		double [][] tmp = new double[4][4];
		for(int i=0;i<4;i++){
			double [] shapeFunctions = setShapeFunctions(this.integralPoints[i].ksi,this.integralPoints[i].eta);
			for(int j=0;j<4;j++){
				tmp[i][j] = shapeFunctions[j];
				System.out.print(tmp[i][j]+"\t");
			}
			System.out.print("\n");
		}
	}

	private double [] setDerivativeShapeFunctions(double val) {
		double [] tmp = new double[4];
		tmp[0] = -(1-val)/4;
		tmp[1] = (1-val)/4 ;
		tmp[2] = (1+val)/4 ;
		tmp[3] = -(1+val)/4 ;
		return tmp;
	}
	//		N1 N2 N3 N4
	//1pc
	//2pc
	//3pc
	//4pc
	double[][] setDNdksi(){ //d shape functions/dksi in integral points
		double [][] tmp = new double[4][4];
		for(int i=0;i<4;i++){
			double [] devShapeFunctions = setDerivativeShapeFunctions(this.integralPoints[i].eta);
			for(int j=0;j<4;j++){
				tmp[i][j] = devShapeFunctions[j];
				System.out.print(tmp[i][j]+"\t");
			}
			System.out.print("\n");
		}
		return tmp;
	}
	double[][] setDNdeta(){ //d shape functions/deta in integral points
		double [][] tmp = new double[4][4];
		for(int i=0;i<4;i++){
			double [] devShapeFunctions = setDerivativeShapeFunctions(this.integralPoints[i].ksi);
			for(int j=0;j<4;j++){
				tmp[i][j] = devShapeFunctions[j];
				System.out.print(tmp[i][j]+"\t");
			}
			System.out.print("\n");
		}
		return tmp;
	}

	public int getNumberOfIntegralPoints() {
		return numberOfIntegralPoints;
	}

	public double[][] getdNdksi() {
		return dNdksi;
	}

	public double[][] getdNdeta() {
		return dNdeta;
	}
}
