//функция: x1^2 + (x3 - 1)^2 + 2*x2 - x1  --> min
//область D: (x1 - 4)^2 + x2^2 + x3^2 <= 4

public class Main {
    static final double EPS1 = 1e-5;
    static final double EPS2 = 1e-10;

    public static void main(String[] args) {
        gradientProjectionMethod();
        System.out.println("\n\n+-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-+" +
                "\n+-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-+\n" +
                "+-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-+\n\n");
        conditionalGradientMethod();
    }

    static void gradientProjectionMethod(){
        Domain d = new Domain();
        Function func = new Function(10, 0.25, 0.001);
        double[] x = new double[]{-2,0.1,3};
        double[] z;

        System.out.println("\t\t\t\t\t\t\t\t\t\tМЕТОД ПРОЕКЦИИ ГРАДИЕНТА\nДанные по каждой итерации:");
        System.out.println("#####   z                                    s           f'(x)                                |x - z|");
        for(int i = 1;;i++) {
            z = findZ(func, d, x);

            while (func.f(x) - func.f(z) < func.sigma * differenceSquare(z, x)) {
                func.s *= func.alpha;
                z = findZ(func, d, x);
            }
            printP(z, x, func, i);

            if (criterionP(z, x, func)) {
                x = z;
            }
            else break;
        }

        System.out.println("________________________________________________________________________________________________________________");
        System.out.println("Результат:");
        System.out.println("     xMin:  [" + x[0] + ", " + x[1] + ", " + x[2] + "]");
        System.out.println("  f(xMin):  " + func.f(x));
    }

    static void conditionalGradientMethod(){
        Domain d = new Domain();
        Function func = new Function(0.451);
        double[] y;
        double[] z = new double[]{4, 0, 0};
        double[] zz = new double[3];
        //z = d.project(z);

        System.out.println("\t\t\t\t\t\t\t\t\t\tМЕТОД УСЛОВНОГО ГРАДИЕНТА\nДанные по каждой итерации:");
        System.out.println("#####   z                                    y                                   " +
                "s          f(zz)      (f'(xk), x - xk)");
        printC(z, new double[]{0,0,0}, func, new double[]{0,0,0}, 0);

        for(int i = 1; ; i++) {
            y = findY(func, d, z);
            changeS(func, z, y);

            zz[0] = z[0] + func.s * (y[0] - z[0]);
            zz[1] = z[1] + func.s * (y[1] - z[1]);
            zz[2] = z[2] + func.s * (y[2] - z[2]);

            printC(z, y, func, zz, i);

            if (criterionC(func, y, z)) {
                z[0] = zz[0];
                z[1] = zz[1];
                z[2] = zz[2];
            } else {
                break;
            }
        }

        System.out.println("________________________________________________________________________________________________________________");
        System.out.println("Результат:");
        System.out.println("     xMin:  [" + zz[0] + ", " + zz[1] + ", " + zz[2] + "]");
        System.out.println("  f(xMin):  " + func.f(zz));
    }

    private static double fk(Function func, double[] yk, double[] zk){
        //fk(x) = (f'(xk), x - xk)
        double res = 0;
        for (int i = 0; i < 3; i++) {
            res += func.dxn(zk, i+1) * (yk[i] - zk[i]);
        }
        return res;
    }

    private static boolean criterionP(double[] z, double[] x, Function func){ //Критерий выхода для метода проекции градиента
        boolean modCriterion = differenceMod(z, x) > EPS1;
        boolean funcCriterion = func.dfMod(x) > EPS1;

        return modCriterion && funcCriterion;
    }
    private static boolean criterionC(Function func, double[] yk, double[] z){ // Критерий выхода для метода условного градиента
        return (fk(func, yk, z) < -EPS1);
    }

    static double differenceMod(double[] z, double[] x){
        double sum = 0;
        for (int i = 0; i < 3; i++) {
            sum += (z[0] - x[0]) * (z[0] - x[0]);
        }
        return Math.sqrt(sum);
    }
    static double differenceSquare(double[] z, double[] x){
        double sum = 0;
        for (int i = 0; i < 3; i++) {
            sum += (z[0] - x[0]) * (z[0] - x[0]);
        }
        return sum;
    }

    static double[] findZ(Function func, Domain d, double[] x){
        double[] newX = new double[3];
        for (int i = 0; i < 3; i++) {
            newX[i] = x[i] - func.s * func.dxn(x, i + 1);
        }
        return d.project(newX);
    }
    static double[] findY(Function func, Domain d, double[] x){
        double[] y = new double[3];
        for (int i = 0; i < 3; i++) {
            y[i] = d.CENTER[i] - ((d.R * func.dxn(x, i + 1))/func.dfMod(x));
        }
        return y;
    }
    private static double phi(Function func, double s, double[] z, double[] y){
        double[] temp = {z[0] + s * (y[0] - z[0]),
                z[1] + s * (y[1] - z[1]),
                z[2] + s * (y[2] - z[2])};
        return func.f(temp);
    }
    static void changeS(Function func, double[] z, double[] y) {
        double a = 0;
        double b = 1;

        double delta = ((b - a) * (3 - Math.sqrt(5))) / 2;
        double x1 = a + delta;
        double x2 = b - delta;

        while (b - a >= 2 * EPS2) {
            if (phi(func, x1, z, y) > phi(func, x2, z, y)){
                a = x1;
                x1 = x2;
                x2 = b + a - x1;
            } else {
                b = x2;
                x2 = x1;
                x1 = a + b - x2;
            }
        }
        func.s = (a + b) / 2;
    }
    //Вспомогательные
    private static void printP(double[] z, double[] x, Function func, int i){
        System.out.printf("%5d   [% 3.6f, % 3.6f, % 3.6f]    %3.6f    [% 3.6f, % 3.6f, % 3.6f]    %3.6f\n", i, z[0], z[1], z[2], func.s, func.dxn(x, 1),
                func.dxn(x, 2), func.dxn(x, 3), differenceMod(z, x));
    }
    private static void printC(double[] z, double[] y, Function func, double[] zz, int i){
        System.out.printf("%5d   [% 3.6f, % 3.6f, % 3.6f]    [% 3.6f, % 3.6f, % 3.6f]   %3.6f   %3.6f   %3.6f\n", i, z[0], z[1], z[2],
                y[0], y[1], y[2], func.s, func.f(zz), fk(func, y, z));
    }
}

final class Domain{
    //Область D: (x1 - 4)^2 + x2^2 + x3^2 <= 4
    //Шар радиуса 2 и с центром в точке (4, 0, 0)
    final double[] CENTER;
    final double R;

    {
        CENTER = new double[]{4.0, 0.0, 0.0};
        R = 2;
    }

    boolean isInDomain(double[] x){
        return (x[0] - 4) * (x[0] - 4) + x[1] * x[1] + x[2] * x[2] <= 4;
    }
    //double checkCircle(double[] x) {return (x[0] - 4) * (x[0] - 4) + x[1] * x[1] + x[2] * x[2];} //если ~= 4, то точка находится на поверхности области

    double[] project(double[] z){
        //функция проекции точки на множество D
        assert z.length == 3 : "Vector length is not equal to 3 / project";
        if(this.isInDomain(z)) return z;

        double[] res = new double[3];
        double zNorm = this.zNorm(z);
        for (int i = 0; i < 3; i++) {
            res[i] = CENTER[i] + R * (z[i] - CENTER[i])/zNorm;
        }

        return res;
    }

    private double zNorm(double[] x){ //z - точка вне области; считает ||z - a||, где a - центр шара
        return Math.sqrt((x[0] - CENTER[0]) * (x[0] - CENTER[0]) +
                (x[1] - CENTER[1]) * (x[1] - CENTER[1]) +
                (x[2] - CENTER[2]) * (x[2] - CENTER[2]));
    }

}

class Function{
    double s;
    double alpha;
    double sigma;

    Function(double s, double alpha, double sigma){
        assert alpha > 0 || alpha < 1 : "Incorrect alpha value";
        assert sigma > 0 || sigma < 1 : "Incorrect sigma value";

        this.s = s;
        this.alpha = alpha;
        this.sigma = sigma;
    }

    Function(double s){
        assert s >= 0 && s <= 1 : "Incorrect s value";
        this.s = s;
    }

    double f(double[] x){
        assert x.length == 3 : "Vector length is not equal to 3 / f";
        return (x[0] * x[0] + (x[2] - 1) * (x[2] - 1) + 2 * x[1] - x[0]);
    }

    private double dx1(double[] x){
        return 2 * x[0] - 1;
    }

    private double dx2(double[] x){
        return 2;
    }

    private double dx3(double[] x){
        return 2 * x[2] - 2;
    }

    double dxn(double[] x, int n){
        switch(n){
            case 1: return dx1(x);
            case 2: return dx2(x);
            case 3: return dx3(x);
            default: throw new IllegalArgumentException();
        }
    }

    double dfMod(double[] x){ // |f'(x)|
        double sum = 0;
        for (int i = 0; i < 3; i++) {
            sum += this.dxn(x, i + 1) * this.dxn(x, i + 1);
        }
        return Math.sqrt(sum);
    }
}


