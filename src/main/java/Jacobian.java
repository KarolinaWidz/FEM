class Jacobian {
	private double upLeft;
	private double upRight;
	private double downLeft;
	private double downRight;


	Jacobian(double dxdksi, double dydksi, double dxdeta, double dydeta) {
		this.upLeft = dxdksi;
		this.upRight = dydksi;
		this.downLeft = dxdeta;
		this.downRight = dydeta;
	}

	void transposedJacobiansDeterminal(){
		this.upLeft = downRight;
		this.upRight = -upRight;
		this.downLeft = -downLeft;
		this.downRight = upLeft;
	}
}
