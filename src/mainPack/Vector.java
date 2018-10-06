package mainPack;

public class Vector extends Matrix{
	Vector(int dimension){
		data = new double[dimension][1];
		rowCount = dimension;
		columnCount = 1;
	}
	
	Vector(double[]arg){
		super(convert(arg));
	}
	
	private static double[][] convert(double[]arg){
		double[][]array = new double[arg.length][1];
		for(int i = 0;i < arg.length;i++) {
			array[i][0] = arg[i];
		}
		return array;
	}
	
	void setValue(double value, int columnIndex) {
		data[0][columnIndex] = value;
	}
	
	double getValue(int columnIndex) {
		return data[0][columnIndex];
	}
	
	double getSize() {
		double sum = 0;
		for(int i = 0;i < data.length;i++) {
			sum += Math.pow(data[i][0], 2);
		}
		return Math.sqrt(sum);
	}
}
