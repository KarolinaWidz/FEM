import static java.lang.Math.sqrt;

class UniversalElement {
	private IntegralPoint [] integralPoints;
	double [] weight;
	double [] shapeFunctions;
	double [][] derivativeShapeFunctionsKsi;
	double [][] derivativeShapeFunctionsEta;


	UniversalElement() {
		this.integralPoints = setIntegralPoints();

	}

	private IntegralPoint[] setIntegralPoints() {
		IntegralPoint [] integralPoints = new IntegralPoint[4];
		integralPoints[0] = new IntegralPoint(-1/sqrt(3),-1/sqrt(3));
		integralPoints[1] = new IntegralPoint(1/sqrt(3),-1/sqrt(3));
		integralPoints[2] = new IntegralPoint(1/sqrt(3),1/sqrt(3));
		integralPoints[3] = new IntegralPoint(-1/sqrt(3),1/sqrt(3));
		return integralPoints;
	}

	private double [] setShapeFunctions(double ksi, double eta) {
		double [] tmp = new double[4];
		tmp[0] = ((1-ksi)*(1-eta))/4;
		tmp[1] = ((1+ksi)*(1-eta))/4 ;
		tmp[2] = ((1+ksi)*(1+eta))/4 ;
		tmp[3] = ((1-ksi)*(1+eta))/4 ;
		return tmp;
	}

	void shapeFunctionInIntegralPoints(){
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

	void dNdksi(){
		double [][] tmp = new double[4][4];
		for(int i=0;i<4;i++){
			double [] devShapeFunctions = setDerivativeShapeFunctions(this.integralPoints[i].eta);
			for(int j=0;j<4;j++){
				tmp[i][j] = devShapeFunctions[j];
				System.out.print(tmp[i][j]+"\t");
			}
			System.out.print("\n");
		}
	}
	void dNdeta(){
		double [][] tmp = new double[4][4];
		for(int i=0;i<4;i++){
			double [] devShapeFunctions = setDerivativeShapeFunctions(this.integralPoints[i].ksi);
			for(int j=0;j<4;j++){
				tmp[i][j] = devShapeFunctions[j];
				System.out.print(tmp[i][j]+"\t");
			}
			System.out.print("\n");
		}

	}

}