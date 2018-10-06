package mainPack;

import java.text.DecimalFormat;
import java.text.NumberFormat;

public class Matrix{
	protected double[][]data;
	protected int rowCount, columnCount;
	public static final int CRAMER_METHOD = 0;
	
	Matrix(){
		
	}
	
	Matrix(int row, int column){
		rowCount = row;
		columnCount = column;
		data = new double[rowCount][columnCount];
	}
	
	Matrix(double[][]arg){
		rowCount = arg.length;
		columnCount = arg[0].length;
		data = new double[rowCount][columnCount];
		setAllValues(arg);
	}
	
	Matrix(Vector[]arg){
		rowCount = arg[0].rowCount;
		columnCount = arg.length;
		data = new double[rowCount][columnCount];
		for(int i = 0;i < columnCount;i++) {
			setColumn(arg[i].getColumn(0),i);
		}
	}
	
	int getRowCount() {
		return rowCount;
	}
	
	int getColumnCount() {
		return columnCount;
	}
	
	void swapRows(int rowIndex1, int rowIndex2) {
		try {
			if(rowIndex1 > -1 && rowIndex1 < rowCount && rowIndex2 > 0 && rowIndex2 < rowCount) {
				double[]swap = data[rowIndex1];
				data[rowIndex1] = data[rowIndex2];
				data[rowIndex2] = swap;
			}else {
				throw new MatrixException("Row couldn't find");
			}
		}catch(MatrixException e) {
			e.printStackTrace();
			System.exit(0);
		}
	}
	
	void swapColumns(int columnIndex1, int columnIndex2) {
		try {
			if(columnIndex1 > -1 && columnIndex1 < columnCount && columnIndex2 > -1 && columnIndex2 < columnCount) {
				for(int i = 0;i < rowCount;i++) {
					double swap = data[i][columnIndex1];
					data[i][columnIndex1] = data[i][columnIndex2];
					data[i][columnIndex2] = swap;
				}
			}else {
				throw new MatrixException("Column couldn't find");
			}
		}catch(MatrixException e) {
			e.printStackTrace();
			System.exit(0);
		}
	}
	
	void setRow(double[]row, int rowIndex) {
		try {
			if(rowIndex > -1 && rowIndex < rowCount && row.length == columnCount) {
				for(int i = 0;i < columnCount;i++) {
					data[rowIndex][i] = row[i];
				}
			}else {
				throw new MatrixException("Row couldn't find or rows doesn't match");
			}
		}catch(MatrixException e) {
			e.printStackTrace();
			System.exit(0);
		}
	}

	double[]getRow(int rowIndex){
		double[]row = null;
		try {
			if(rowIndex > -1 && rowIndex < rowCount) {
				row = new double[columnCount];
				for(int i = 0;i < columnCount;i++) {
					row[i] = data[rowIndex][i];
				}
			}else {
				throw new MatrixException("Row couldn't find");
			}
		}catch(MatrixException e) {
			e.printStackTrace();
			System.exit(0);
		}
		return row;
	}
	
	void setColumn(double[]column, int columnIndex) {
		try {
			if(columnIndex > -1 && columnIndex < columnCount && column.length == rowCount) {
				for(int i = 0;i < rowCount;i++) {
					data[i][columnIndex] = column[i];
				}
			}else {
				throw new MatrixException("Columm couldn't find or columns doesn't match");
			}
		}catch(MatrixException e) {
			e.printStackTrace();
			System.exit(0);
		}
	}
	
	double[]getColumn(int columnIndex){
		double[]column = null;
		try {
			if(columnIndex > -1 && columnIndex < columnCount) {
				column = new double[rowCount];
				for(int i = 0;i < rowCount;i++) {
					column[i] = data[i][columnIndex];
				}
			}else {
				throw new MatrixException("Column couldn't find");
			}
		}catch(MatrixException e) {
			e.printStackTrace();
			System.exit(0);
		}
		return column;
	}
	
	void setValue(double value, int rowIndex, int columnIndex) {
		data[rowIndex][columnIndex] = value;
	}
	
	double getValue(int rowIndex, int columnIndex) {
		return data[rowIndex][columnIndex];
	}
	
	Matrix multiplicate(double n) {
		Matrix m1 = new Matrix(rowCount, columnCount);
		for(int i = 0;i < rowCount;i++) {
			for(int j = 0;j < columnCount;j++) {
				m1.setValue(n*getValue(i,j), i, j);
			}
		}
		return m1;
	}
	
	void setAllValues(double[][]matrix) {
		try {
			if(isSameDimension(matrix)) {
				for(int i = 0;i < rowCount;i++) {
					for(int j = 0;j < columnCount;j++) {
						setValue(matrix[i][j], i, j);
					}
				}
			}else {
				throw new MatrixException("Matrix dimension isn't equal");
			}
		}catch(MatrixException e) {
			e.printStackTrace();
			System.exit(0);
		}
	}
	
	double[][] getAllValues() {
		return data;
	}
	
	boolean isSameDimension(Matrix matrix) {
		if(matrix.getRowCount() == getRowCount() && matrix.getColumnCount() == getColumnCount()) {
			return true;
		}else {
			return false;
		}
	}
	
	boolean isSameDimension(double[][] matrix) {
		if(matrix.length == getRowCount() && matrix[0].length == getColumnCount()) {
			return true;
		}else {
			return false;
		}
	}
	
	static Matrix sum(Matrix m1, Matrix m2) {
		try {
			if(m1.isSameDimension(m2)) {
				Matrix m3 = new Matrix(m1.getRowCount(),m1.getColumnCount());
				for(int i = 0;i < m1.getRowCount();i++) {
					for(int j = 0;j < m1.getColumnCount();j++) {
						m3.setValue(m1.getValue(i, j) + m2.getValue(i, j), i, j);
					}
				}
				return m3;
			}else {
				throw new MatrixException();
			}
		}catch(MatrixException e) {
			e.printStackTrace();
			System.exit(0);
		}
		return null;
	}
	
	static Matrix ex(Matrix m1, Matrix m2) {
		try {
			if(m1.isSameDimension(m2)) {
				Matrix m3 = new Matrix(m1.getRowCount(),m1.getColumnCount());
				for(int i = 0;i < m1.getRowCount();i++) {
					for(int j = 0;j < m1.getColumnCount();j++) {
						m3.setValue(m1.getValue(i, j) - m2.getValue(i, j), i, j);
					}
				}
				return m3;
			}else {
				throw new MatrixException();
			}
		}catch(MatrixException e) {
			e.printStackTrace();
			System.exit(0);
		}
		return null;
	}
	
	static Matrix multiplicate(Matrix m1, Matrix m2) {
		try {
			if(m1.getColumnCount() == m2.getRowCount()) {
				double sum = 0;
				Matrix m3 = new Matrix(m1.getRowCount(),m2.getColumnCount());
				for(int i = 0;i < m3.getRowCount();i++) {
					for(int j = 0;j < m3.getColumnCount();j++) {
						for(int w = 0;w < m1.getColumnCount();w++) {
							sum += m1.getValue(i, w) * m2.getValue(w, j);
						}
						m3.setValue(sum, i, j);
						sum = 0;
					}
				}
				return m3;
			}else {
				throw new MatrixException();
			}
		}catch(MatrixException e) {
			e.printStackTrace();
			System.exit(0);
		}
		return null;
	}
	
	public double determinant() {
		try {
			if(rowCount != columnCount) {
				throw new MatrixException("It should be square matrix");
			}else {
				if(rowCount == 1) {
					return getValue(0,0);
				}else {
					double sum = 0;
					boolean sign = true;
					for(int i = 0;i < columnCount;i++) {
						sum += (sign)? getValue(0,i)*this.minor(0, i) : -getValue(0,i)*this.minor(0, i);
						sign = !sign;
					}
					return sum;
				}
			}
		}catch(MatrixException e) {
			e.printStackTrace();
			System.exit(0);
		}
		return 0;
	}
	
	public double minor(int rowIndex, int columnIndex) {
		try {
			if(rowIndex > -1 && rowIndex < rowCount && columnIndex > -1 && columnIndex < columnCount) {
				Matrix matrix = new Matrix(rowCount - 1,columnCount - 1);
				int row = 0,column = 0;
				for(int i = 0;i < rowCount;i++) {
					if(i == rowIndex)continue;
					for(int j = 0;j < columnCount;j++) {
						if(j == columnIndex)continue;
						matrix.setValue(getValue(i,j), row, column);
						column++;
					}
					column = 0;
					row++;
				}
				return matrix.determinant();
			}else {
				throw new MatrixException("Out of row or column");
			}
		}catch(MatrixException e) {
			e.printStackTrace();
			System.exit(0);
		}
		return 0;
	}
	
	public Matrix cofactor() {
		Matrix matrix = minorMatrix();
		for(int i = 0;i < rowCount;i++) {
			for(int j = 0;j < columnCount;j++) {
				matrix.setValue(Math.pow(-1, i + j)*matrix.getValue(i, j), i, j);
			}
		}
		return matrix;
	}
	
	public Matrix transpose() {
		Matrix matrix = new Matrix(columnCount,rowCount);
		for(int i = 0;i < rowCount;i++) {
			for(int j = 0;j < columnCount;j++) {
				matrix.setValue(getValue(i,j), j, i);
			}
		}
		return matrix;
	}
	
	public Matrix minorMatrix() {
		Matrix matrix = new Matrix(rowCount,columnCount);
		for(int i = 0;i < rowCount;i++) {
			for(int j = 0;j < columnCount;j++) {
				matrix.setValue(minor(i,j), i, j);
			}
		}
		return matrix;
	}
	
	public Matrix inverse() {
		if(rowCount == 1) {
			return new Matrix(new double[][] {{1.0 / this.getValue(0, 0)}});
		}
		return cofactor().transpose().multiplicate(1 / determinant());
	}
	
	public static Matrix generateUnitMatrix(int dimension) {
		Matrix matrix = new Matrix(dimension,dimension);
		for(int i = 0;i < dimension;i++) {
			for(int j = 0;j < dimension;j++) {
				matrix.setValue((i == j)? 1 : 0, i, j);
			}
		}
		return matrix;
	}
	
	public static Matrix solve(Matrix matrix1, Matrix matrix2) {
		if(matrix1.determinant() == 0) {
			System.out.println("infinite solutions");
			System.exit(0);
		}
		return Matrix.multiplicate(matrix1.inverse(),matrix2);
	}
	
	public Matrix getGramSchmidt() {
		Vector[]bases = new Vector[columnCount];
		for(int i = 0;i < columnCount;i++) {
			bases[i] = new Vector(getColumn(i));
			for(int k = 0;k < i;k++) {
				Vector v1 = new Vector(getColumn(k));
				bases[i] = new Vector(Matrix.ex(bases[i], Matrix.multiplicate(v1.getProjectionMatrix(), bases[i])).getColumn(0));
			}
		}
		return new Matrix(bases);
	}
	
	public Matrix getProjectionMatrix() {
		return Matrix.multiplicate(this, Matrix.multiplicate(Matrix.multiplicate(this.transpose(), this).inverse(), this.transpose()));
	}
	public Matrix getLeastSquareMatrix() {
		return Matrix.multiplicate(Matrix.multiplicate(this.transpose(), this).inverse(), this.transpose());
	}
	
	public boolean isEqual(Matrix matrix) {
		if(rowCount == matrix.getRowCount() && columnCount == matrix.getColumnCount()) {
			for(int i = 0;i < rowCount;i++) {
				for(int j = 0;j < columnCount;j++) {
					if(getValue(i, j) != matrix.getValue(i, j))return false;
				}
			}
			return true;
		}else {
			return false;
		}
	}
	
	@Override
	public String toString() {
		StringBuffer sf = new StringBuffer();
		NumberFormat formatter = new DecimalFormat("#.###");   
		for(int i = 0;i < rowCount;i++) {
			sf.append("[");
			for(int j = 0;j < columnCount;j++) {
				sf.append(formatter.format(getValue(i,j)));
				if(j < columnCount - 1) {
					sf.append(", ");
				}
			}
			sf.append("]" + System.lineSeparator());
		}
		return sf.toString();
	}
}