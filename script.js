// Example Nonlinear Finite Difference Method with full detailed outputs

// Define f(x, y, y') and its partial derivatives fy, fyp
function f(x, y, yp) {
    return y; // Example: f = y
}

function fy(x, y, yp) {
    return 1; // df/dy = 1
}

function fyp(x, y, yp) {
    return 0; // df/dy' = 0
}

// Parameters
const a = 0, b = 1;
const alpha = 0, beta = 0;
const N = 4;
const TOL = 1e-5;
const M = 10;

// Initialize
const h = (b - a) / (N + 1);
let w = new Array(N + 2);

// Set boundary conditions
w[0] = alpha;
w[N+1] = beta;

// Initial guess (linear interpolation)
for (let i = 1; i <= N; i++) {
    w[i] = alpha + (i * (beta - alpha) * h) / (b - a);
}

// Show initial values
console.log(`h = ${h}`);
console.log("Initial guess:");
for (let i = 0; i <= N+1; i++) {
    console.log(`w[${i}] = ${w[i]}`);
}

// Start iterations
let k = 1;
while (k <= M) {
    console.log(`\nIteration ${k}:`);
    let aArr = [], bArr = [], cArr = [], dArr = [];

    // Step 5
    let x = a + h;
    let t = (w[2] - alpha) / (2 * h);
    let a1 = 2 + h**2 * fy(x, w[1], t);
    let b1 = -1 + (h/2) * fyp(x, w[1], t);
    let d1 = -(2*w[1] - w[2] - alpha + h**2*f(x, w[1], t));
    aArr[1] = a1; bArr[1] = b1; dArr[1] = d1;

    console.log(`i=1, x=${x.toFixed(4)}, t=${t.toFixed(4)}, a1=${a1.toFixed(4)}, b1=${b1.toFixed(4)}, d1=${d1.toFixed(4)}`);

    // Step 6
    for (let i = 2; i <= N-1; i++) {
        x = a + i * h;
        t = (w[i+1] - w[i-1]) / (2 * h);
        let ai = 2 + h**2 * fy(x, w[i], t);
        let bi = -1 + (h/2) * fyp(x, w[i], t);
        let ci = -1 - (h/2) * fyp(x, w[i], t);
        let di = -(2*w[i] - w[i+1] - w[i-1] + h**2*f(x, w[i], t));
        aArr[i] = ai; bArr[i] = bi; cArr[i] = ci; dArr[i] = di;

        console.log(`i=${i}, x=${x.toFixed(4)}, t=${t.toFixed(4)}, ai=${ai.toFixed(4)}, bi=${bi.toFixed(4)}, ci=${ci.toFixed(4)}, di=${di.toFixed(4)}`);
    }

    // Step 7
    x = b - h;
    t = (beta - w[N-1]) / (2*h);
    let aN = 2 + h**2 * fy(x, w[N], t);
    let cN = -1 - (h/2) * fyp(x, w[N], t);
    let dN = -(2*w[N] - w[N-1] - beta + h**2*f(x, w[N], t));
    aArr[N] = aN; cArr[N] = cN; dArr[N] = dN;

    console.log(`i=${N}, x=${x.toFixed(4)}, t=${t.toFixed(4)}, aN=${aN.toFixed(4)}, cN=${cN.toFixed(4)}, dN=${dN.toFixed(4)}`);

    // Step 8-10: Solve the tridiagonal system
    let l = [], u = [], z = [];
    l[1] = aArr[1];
    u[1] = bArr[1] / l[1];
    z[1] = dArr[1] / l[1];
    console.log(`Forward Elimination i=1: l1=${l[1].toFixed(4)}, u1=${u[1].toFixed(4)}, z1=${z[1].toFixed(4)}`);

    for (let i = 2; i <= N-1; i++) {
        l[i] = aArr[i] - cArr[i]*u[i-1];
        u[i] = bArr[i]/l[i];
        z[i] = (dArr[i] - cArr[i]*z[i-1])/l[i];
        console.log(`Forward Elimination i=${i}: li=${l[i].toFixed(4)}, ui=${u[i].toFixed(4)}, zi=${z[i].toFixed(4)}`);
    }

    l[N] = aArr[N] - cArr[N]*u[N-1];
    z[N] = (dArr[N] - cArr[N]*z[N-1])/l[N];
    console.log(`Forward Elimination i=${N}: lN=${l[N].toFixed(4)}, zN=${z[N].toFixed(4)}`);

    // Step 11: Back substitution
    let v = [];
    v[N] = z[N];
    w[N] += v[N];
    console.log(`Back Substitution i=${N}: vN=${v[N].toFixed(4)}, wN=${w[N].toFixed(4)}`);

    for (let i = N-1; i >= 1; i--) {
        v[i] = z[i] - u[i]*v[i+1];
        w[i] += v[i];
        console.log(`Back Substitution i=${i}: vi=${v[i].toFixed(4)}, wi=${w[i].toFixed(4)}`);
    }

    // Step 13: Check convergence
    let maxV = Math.max(...v.slice(1).map(Math.abs));
    console.log(`Max correction = ${maxV.toExponential(3)}`);
    if (maxV < TOL) {
        console.log("\nConvergence achieved.\n");
        break;
    }

    k++;
}

if (k > M) {
    console.log("\nMaximum number of iterations exceeded.");
}

// Step 14: Output final table
console.log("\nFinal Table:");
for (let i = 0; i <= N+1; i++) {
    let xi = a + i*h;
    console.log(`x=${xi.toFixed(4)}, w=${w[i].toFixed(6)}`);
}

// Adjust input width as user types
function adjustWidth(input) {
    input.style.width = ((input.value.length + 1) * 8) + 'px';
}
