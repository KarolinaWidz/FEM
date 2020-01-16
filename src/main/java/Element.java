import org.ejml.simple.SimpleMatrix;
import java.util.Arrays;
import java.lang.*;
import static java.lang.Math.pow;
import static java.lang.Math.sqrt;

/**
 * Class with global finite element
 */
class Element {
	private int elementID;
	private int [] nodesID;
	private Node [] nodes;
	private UniversalElement universalElement;
	private SimpleMatrix hMatrix;
	private SimpleMatrix cMatrix;
	private double [] jacobianDet;
	private SimpleMatrix hbcMatrix;
	private Jacobian [] jacobian;
	private Jacobian [] cofactorJacobian; //macierz dopelnien algebraicznych

	Element(int elementID, int nodesPerHeight) {
		this.elementID = elementID;
		this.nodesID = setNodesID(elementID,nodesPerHeight);
		this.universalElement = new UniversalElement();
		this.jacobian = new Jacobian[universalElement.getNumberOfIntegralPoints()];
		this.jacobianDet = new double[universalElement.getNumberOfIntegralPoints()];
		this.cofactorJacobian = new Jacobian[universalElement.getNumberOfIntegralPoints()];
		this.hMatrix = new SimpleMatrix(universalElement.getNumberOfIntegralPoints(),universalElement.getNumberOfIntegralPoints());
		this.cMatrix = new SimpleMatrix(universalElement.getNumberOfIntegralPoints(),universalElement.getNumberOfIntegralPoints());
		this.hbcMatrix = new SimpleMatrix(universalElement.getNumberOfIntegralPoints(),universalElement.getNumberOfIntegralPoints());
	}

	private int [] setNodesID(int elementID, int nodesPerHeight){
		nodesID = new int [4];
		if(elementID<5) {
			nodesID[0] = elementID;
			nodesID[1] = nodesID[0] + nodesPerHeight;
			nodesID[2] = nodesID[1] + 1;
			nodesID[3] = nodesID[0] + 1;
		}
		else if(elementID<10) {
			nodesID[0] = elementID + 1;
			nodesID[1] = nodesID[0] + nodesPerHeight;
			nodesID[2] = nodesID[1] + 1;
			nodesID[3] = nodesID[0] + 1;
		}
		else {
			nodesID[0] = elementID + 2;
			nodesID[1] = nodesID[0] + nodesPerHeight;
			nodesID[2] = nodesID[1] + 1;
			nodesID[3] = nodesID[0] + 1;
		}
		return nodesID;
	}

	void calculateJacobianDeterminal(){
		double [] dxdksi = new double[4] ;
		double [] dydksi = new double[4] ;
		double [] dxdeta = new double[4] ;
		double [] dydeta = new double[4] ;
		double [][] dNdksi = universalElement.getdNdksi();
		double [][] dNdeta = universalElement.getdNdeta();

		for(int i=0;i<universalElement.getNumberOfIntegralPoints();i++){
			dxdksi[i] = calculateDerivativeX(dNdksi[i]);
			dydksi[i] = calculateDerivativeY(dNdksi[i]);
			dxdeta[i] = calculateDerivativeX(dNdeta[i]);
			dydeta[i] = calculateDerivativeY(dNdeta[i]);
			this.jacobian[i] = new Jacobian(dxdksi[i], dydksi[i], dxdeta[i], dydeta[i]);
			this.jacobianDet[i] = dxdksi[i]*dydeta[i]-dydksi[i]*dxdeta[i];
			this.cofactorJacobian[i] = jacobian[i].transposedJacobiansDeterminal();

		}
	}


	/**
	 * Method used in calculating jacobians
	 * derivative - dNdksi or dNdeta
	 */
	private double calculateDerivativeX(double[] derivative){
		double tmp;
			tmp = derivative[0]*nodes[0].getX()+
					derivative[1]*nodes[1].getX()+
					derivative[2]*nodes[2].getX()+
					derivative[3]*nodes[3].getX();

		return tmp;
	}

	/**
	 * Method used in calculating jacobians
	 * derivative - dNdksi or dNdeta
	 */
	private double calculateDerivativeY(double[] derivative){
		double tmp;
		tmp = derivative[0]*nodes[0].getY()+
				derivative[1]*nodes[1].getY()+
				derivative[2]*nodes[2].getY()+
				derivative[3]*nodes[3].getY();

		return tmp;
	}

	private double dNdx(int integralPoint, int numberOfShapeFunction){

		return (1/this.jacobianDet[integralPoint])*((this.cofactorJacobian[integralPoint].getUpLeft()
				*universalElement.getdNdksi()[integralPoint][numberOfShapeFunction])-(this.cofactorJacobian[integralPoint].getUpRight()
				*universalElement.getdNdeta()[integralPoint][numberOfShapeFunction]));
	}

	private double dNdy(int integralPoint, int numberOfShapeFunction){
		return (1/ this.jacobianDet[integralPoint])*((this.cofactorJacobian
				[integralPoint].getDownRight()*universalElement.getdNdeta()[integralPoint][numberOfShapeFunction
				])-(this.cofactorJacobian[integralPoint].getDownLeft()
				*universalElement.getdNdksi()[integralPoint][numberOfShapeFunction]));
	}

	void calculateHLocalMatrix(double conductivity) {
		SimpleMatrix [] hLocal = new SimpleMatrix[universalElement.getNumberOfIntegralPoints()];
		for (int i = 0; i < universalElement.getNumberOfIntegralPoints(); i++) {
			SimpleMatrix dNdx = dvectorNdx(i);
			SimpleMatrix dNdy = dvectorNdy(i);
			SimpleMatrix dx = dNdx.mult(dNdx.transpose());
			SimpleMatrix dy = dNdy.mult(dNdy.transpose());
			dx = dx.scale(this.jacobianDet[i]);
			dy = dy.scale(this.jacobianDet[i]);
			hLocal[i] = dx.plus(dy);
		}

		for (int i = 0; i < universalElement.getNumberOfIntegralPoints(); i++) {
			hLocal[i] = hLocal[i].scale(conductivity*universalElement.getIntegralPoints()[i].
					ksiWeight*universalElement.getIntegralPoints()[i].etaWeight);
			this.hMatrix = this.hMatrix.plus(hLocal[i]);
		}
		System.out.println(this.hMatrix);

	}

	private SimpleMatrix dvectorNdx(int integralPoint){
		double [][] tmp = new double[universalElement.getNumberOfIntegralPoints()][1];
		for(int i=0;i<universalElement.getNumberOfIntegralPoints();i++){
			tmp[i][0]=dNdx(integralPoint,i);
		}
		return new SimpleMatrix(tmp);
	}

	private SimpleMatrix dvectorNdy(int integralPoint){
		double [][] tmp = new double[universalElement.getNumberOfIntegralPoints()][1];
		for(int i=0;i<universalElement.getNumberOfIntegralPoints();i++){
			tmp[i][0]=dNdy(integralPoint,i);
		}
		return new SimpleMatrix(tmp);

	}


	/**
	 * specificHeat - c
	 * density - ro
	 */
	void calculateCLocalMatrix(int specificHeat, int density){
		SimpleMatrix [] cLocal = new SimpleMatrix[universalElement.getNumberOfIntegralPoints()];
		for(int i=0;i<universalElement.getNumberOfIntegralPoints();i++){
			double [] tmp = universalElement.shapeFunctions(universalElement.getIntegralPoints()[i].ksi,universalElement.getIntegralPoints()[i].eta);
			double[][] tmpToSimpleMatrix = new double[universalElement.getNumberOfIntegralPoints()][1];
			for(int j=0;j<universalElement.getNumberOfIntegralPoints();j++){
				tmpToSimpleMatrix[j][0] = tmp[j];
			}
			SimpleMatrix shapeFunctions = new SimpleMatrix(tmpToSimpleMatrix);
			SimpleMatrix shapeFunctionsTransposed = shapeFunctions.transpose();
			SimpleMatrix shapeFunctionsMultiplication = shapeFunctions.mult(shapeFunctionsTransposed);
			double scalar = density * specificHeat * this.jacobianDet[i];
			cLocal[i] = shapeFunctionsMultiplication.scale(scalar);
		}

		for(int i=0;i<universalElement.getNumberOfIntegralPoints();i++){
			cLocal[i] = cLocal[i].scale(universalElement.getIntegralPoints()[i].ksiWeight*universalElement.getIntegralPoints()[i].etaWeight);
			this.cMatrix = this.cMatrix.plus(cLocal[i]);
		}
		System.out.println(this.cMatrix);
	}

	void calculateHBC(int alfa){

		double [] length = new double [4];
		double [] x = new double[4];
		double [] y = new double[4];
		for(int i=0;i<universalElement.getNumberOfIntegralPoints();i++){
			x[i] = getNodes()[i].getX();
			y[i] = getNodes()[i].getY();
		}
		length[0] = sqrt(pow((x[1]-x[0]),2)+pow((y[1]-y[0]),2));
		length[1] = sqrt(pow((x[1]-x[2]),2)+pow((y[1]-y[2]),2));
		length[2] = sqrt(pow((x[2]-x[3]),2)+pow((y[2]-y[3]),2));
		length[3] = sqrt(pow((x[3]-x[0]),2)+pow((y[3]-y[0]),2));

		Surface surface [] = new Surface[4];
		surface[0] = new Surface(-1/sqrt(3),-1,1/sqrt(3),-1,length[0]/2);
		surface[1] = new Surface(1,-1/sqrt(3),1,1/sqrt(3),length[1]/2);
		surface[2] = new Surface(1/sqrt(3),1,-1/sqrt(3),1,length[2]/2);
		surface[3] = new Surface(-1,1/sqrt(3),-1,-1/sqrt(3),length[3]/2);

		double [][][] tmpFirstPoint = new double[4][1][4];
		double [][][] tmpSecondPoint = new double[4][1][4];
		SimpleMatrix[] shapeFunctionsInFirstPoint = new SimpleMatrix[4];
		SimpleMatrix[] shapeFunctionsInSecondPoint = new SimpleMatrix[4];
		SimpleMatrix[] resultInFirstPoint = new SimpleMatrix[4];
		SimpleMatrix[] resultInSecondPoint = new SimpleMatrix[4];
		SimpleMatrix sum [] = new SimpleMatrix[4];
		for(int surfaceNum=0;surfaceNum<4;surfaceNum++){
			for(int pointNum=0;pointNum<2;pointNum++){
				for(int functionNum=0;functionNum<4;functionNum++){
					tmpFirstPoint[surfaceNum][0][functionNum]=universalElement.shapeFunctions
							(surface[surfaceNum].getPointsOnSurface()[0].ksi,surface
									[surfaceNum].getPointsOnSurface()[0].eta)[functionNum];
					tmpSecondPoint[surfaceNum][0][functionNum]=universalElement.shapeFunctions
							(surface[surfaceNum].getPointsOnSurface()[1].ksi,surface
									[surfaceNum].getPointsOnSurface()[1].eta)[functionNum];
				}
			}
		}
		for(int surfaceNum=0;surfaceNum<4;surfaceNum++){
			shapeFunctionsInFirstPoint[surfaceNum] = new SimpleMatrix(tmpFirstPoint[surfaceNum]);
			shapeFunctionsInSecondPoint[surfaceNum] = new SimpleMatrix(tmpSecondPoint[surfaceNum]);
			resultInFirstPoint[surfaceNum] = shapeFunctionsInFirstPoint[surfaceNum].transpose().mult(shapeFunctionsInFirstPoint[surfaceNum]).scale(alfa);
			resultInSecondPoint[surfaceNum] = shapeFunctionsInSecondPoint[surfaceNum].transpose().mult(shapeFunctionsInSecondPoint[surfaceNum]).scale(alfa);
			sum[surfaceNum] = resultInFirstPoint[surfaceNum].plus(resultInSecondPoint[surfaceNum]).scale(surface[surfaceNum].detJ);

		}

		SimpleMatrix tmpResult1;
		SimpleMatrix tmpResult2;
		SimpleMatrix tmpResult3;
		SimpleMatrix tmpResult4;
		//tmpResult1 = sum[0].scale(getNodes()[0].getIntBoundaryCondition()*getNodes()[1].getIntBoundaryCondition());
		//tmpResult2 = sum[1].scale(getNodes()[1].getIntBoundaryCondition()*getNodes()[2].getIntBoundaryCondition());
		//tmpResult3 = sum[2].scale(getNodes()[2].getIntBoundaryCondition()*getNodes()[3].getIntBoundaryCondition());
		//tmpResult4 = sum[3].scale(getNodes()[3].getIntBoundaryCondition()*getNodes()[0].getIntBoundaryCondition());

		//this.hbcMatrix = tmpResult1.plus(tmpResult2).plus(tmpResult3).plus(tmpResult4);
		//System.out.println(this.hbcMatrix);

	}


	public int getElementID() {
		return elementID;
	}

	public void setElementID(int elementID) {
		this.elementID = elementID;
	}

	int[] getNodesID() {
		return nodesID;
	}

	public void setNodesID(int[] nodesID) {
		this.nodesID = nodesID;
	}

	void setNodes(Node[] nodes) {
		this.nodes = nodes;
	}

	public Node[] getNodes() {
		return nodes;
	}

	public int[] printAllNodes() {
		int [] tmp  = new int[4];
		for(int i=0;i<4;i++){
			tmp[i]=nodes[i].getNodeID();
		}
		return tmp;
}

	public String toStringNode(){
		return Arrays.toString(nodesID);
	}

	@Override
	public String toString() {
		return "Element{" +
				"elementID=" + elementID +
				", nodesID=" + Arrays.toString(nodesID) +
				'}';
	}
}
