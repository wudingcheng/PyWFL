void add(double a, double b, double *c) {
    *c = a + b;
}

int main() {
    double a, b, c;
    a = 1.0;
    b = 2.0;
    add(a, b, &c);
    printf("%f + %f = %f\n", a, b, c);
}