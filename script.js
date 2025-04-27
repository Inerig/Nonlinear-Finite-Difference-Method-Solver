function parseFunction(funcStr) {
    try {
        // Validate the function string
        if (!/^[a-zA-Z0-9\s\+\-\*\/\(\)\.\,\[\]\^\%\&\|\!\<\>\=\~\`\@\#\$\^\\]+$/.test(funcStr)) {
            throw new Error("Invalid characters in function.");
       function adjustWidth(input) {
    input.style.width = (input.value.length + 2) + 'ch';
}

function fyApprox(f, x, y, yPrime) {
    const delta = 1e-8;
    const f1 = f(x, y + delta, yPrime);
    const f2 = f(x, y - delta, yPrime);
    return (f1 - f2) / (2 * delta); // Approximate partial derivative fy
}

function solveAndDisplay() {
    const fInput = document.getElementById('functionInput').value;
    const a = parseFloat(document.getElementById('a').value);
    const b = parseFloat(document.getElementById('b').value);
    const alpha = parseFloat(document.getElementById('alpha').value);
    const beta = parseFloat(document.getElementById('beta').value);
    const N = parseInt(document.getElementById('n').value);
    const TOL = parseFloat(document.getElementById('tolerance').value);
    const M = parseInt(document.getElementById('maxIterations').value);
    const decimals = parseInt(document.getElementById('decimalPlaces').value);

    const h = (b - a) / (N + 1);
    const x = [];
    for (let i = 0; i <= N + 1; i++) {
        x.push(a + i * h);
    }

    const f = new Function('x', 'y', 'yPrime', `return ${fInput};`);

    let w = [];
    w[0] = alpha;
    for (let i = 1; i <= N; i++) {
        w[i] = alpha + (i * (beta - alpha) / (b - a)) * h;
    }
    w[N + 1] = beta;

    const calculationsDiv = document.getElementById('nldfCalculationsContent');
    const tbody = document.getElementById('nldfResultTable').querySelector('tbody');
    calculationsDiv.innerHTML = '';
    tbody.innerHTML = '';

    let success = false;
    let k = 1;
    while (k <= M) {
        const aCoeff = [], bCoeff = [], cCoeff = [], dCoeff = [];
        calculationsDiv.innerHTML += `<div class="step"><strong>Iteration ${k}:</strong></div>`;

        // Step 5
        let xi = a + h;
        let t = (w[2] - alpha) / (2 * h);
        aCoeff[1] = 2 + h * h * fyApprox(f, xi, w[1], t);
        bCoeff[1] = -1 + (h / 2) * fyApprox(f, xi, w[1], t);
        dCoeff[1] = -(2 * w[1] - w[2] - alpha + h * h * f(xi, w[1], t));

        calculationsDiv.innerHTML += `<div>Step 5: x = ${xi.toFixed(decimals)}, t = ${t.toFixed(decimals)}, a1 = ${aCoeff[1].toFixed(decimals)}, b1 = ${bCoeff[1].toFixed(decimals)}, d1 = ${dCoeff[1].toFixed(decimals)}</div>`;

        // Step 6
        for (let i = 2; i <= N - 1; i++) {
            xi = a + i * h;
            t = (w[i + 1] - w[i - 1]) / (2 * h);
            aCoeff[i] = 2 + h * h * fyApprox(f, xi, w[i], t);
            bCoeff[i] = -1 + (h / 2) * fyApprox(f, xi, w[i], t);
            cCoeff[i] = -1 - (h / 2) * fyApprox(f, xi, w[i], t);
            dCoeff[i] = -(2 * w[i] - w[i + 1] - w[i - 1] + h * h * f(xi, w[i], t));

            calculationsDiv.innerHTML += `<div>Step 6: i = ${i}, x = ${xi.toFixed(decimals)}, t = ${t.toFixed(decimals)}, ai = ${aCoeff[i].toFixed(decimals)}, bi = ${bCoeff[i].toFixed(decimals)}, ci = ${cCoeff[i].toFixed(decimals)}, di = ${dCoeff[i].toFixed(decimals)}</div>`;
        }

        // Step 7
        xi = b - h;
        t = (beta - w[N - 1]) / (2 * h);
        aCoeff[N] = 2 + h * h * fyApprox(f, xi, w[N], t);
        cCoeff[N] = -1 - (h / 2) * fyApprox(f, xi, w[N], t);
        dCoeff[N] = -(2 * w[N] - w[N - 1] - beta + h * h * f(xi, w[N], t));

        calculationsDiv.innerHTML += `<div>Step 7: x = ${xi.toFixed(decimals)}, t = ${t.toFixed(decimals)}, aN = ${aCoeff[N].toFixed(decimals)}, cN = ${cCoeff[N].toFixed(decimals)}, dN = ${dCoeff[N].toFixed(decimals)}</div>`;

        // Step 8 to 12 (Thomas algorithm for tridiagonal system)
        const l = [], u = [], z = [], v = [];

        l[1] = aCoeff[1];
        u[1] = bCoeff[1] / l[1];
        z[1] = dCoeff[1] / l[1];

        for (let i = 2; i <= N - 1; i++) {
            l[i] = aCoefunction solveAndDisplay() {
    // Read input values
    const funcInput = document.getElementById('functionInput').value;
    const f = new Function('x', 'y', 'yPrime', `return ${funcInput};`);
    const a = parseFloat(document.getElementById('a').value);
    const b = parseFloat(document.getElementById('b').value);
    const alpha = parseFloat(document.getElementById('alpha').value);
    const beta = parseFloat(document.getElementById('beta').value);
    const N = parseInt(document.getElementById('n').value);
    const TOL = parseFloat(document.getElementById('tolerance').value);
    const M = parseInt(document.getElementById('maxIterations').value);
    const decimalPlaces = parseInt(document.getElementById('decimalPlaces').value);

    const h = (b - a) / (N + 1);
    const w = new Array(N + 2);
    w[0] = alpha;
    w[N + 1] = beta;

    // Initial guess for w1, w2, ..., wN (linear interpolation)
    for (let i = 1; i <= N; i++) {
        w[i] = alpha + (i * (beta - alpha) * h) / (b - a);
    }

    let k = 1;
    const calculationsDiv = document.getElementById('nldfCalculationsContent');
    calculationsDiv.innerHTML = ""; // clear previous calculations
    const resultTable = document.getElementById('nldfResultTable').querySelector('tbody');
    resultTable.innerHTML = "";

    while (k <= M) {
        const aArray = new Array(N + 1);
        const bArray = new Array(N + 1);
        const cArray = new Array(N + 1);
        const dArray = new Array(N + 1);

        // Step 5
        let x = a + h;
        let t = (w[2] - alpha) / (2 * h);
        aArray[1] = 2 + Math.pow(h, 2) * fy(f, x, w[1], t);
        bArray[1] = -1 + (h / 2) * fy(f, x, w[1], t);
        dArray[1] = -(2 * w[1] - w[2] - alpha + Math.pow(h, 2) * f(x, w[1], t));

        // Step 6
        for (let i = 2; i <= N - 1; i++) {
            x = a + i * h;
            t = (w[i + 1] - w[i - 1]) / (2 * h);
            aArray[i] = 2 + Math.pow(h, 2) * fy(f, x, w[i], t);
            bArray[i] = -1 + (h / 2) * fy(f, x, w[i], t);
            cArray[i] = -1 - (h / 2) * fy(f, x, w[i], t);
            dArray[i] = -(2 * w[i] - w[i + 1] - w[i - 1] + Math.pow(h, 2) * f(x, w[i], t));
        }

        // Step 7
        x = b - h;
        t = (beta - w[N - 1]) / (2 * h);
        aArray[N] = 2 + Math.pow(h, 2) * fy(f, x, w[N], t);
        cArray[N] = -1 - (h / 2) * fy(f, x, w[N], t);
        dArray[N] = -(2 * w[N] - w[N - 1] - beta + Math.pow(h, 2) * f(x, w[N], t));

        // Step 8
        const l = new Array(N + 1);
        const u = new Array(N + 1);
        const z = new Array(N + 1);

        l[1] = aArray[1];
        u[1] = bArray[1] / l[1];
        z[1] = dArray[1] / l[1];

        // Step 9
        for (let i = 2; i <= N - 1; i++) {
            l[i] = aArray[i] - cArray[i] * u[i - 1];
            u[i] = bArray[i] / l[i];
            z[i] = (dArray[i] - cArray[i] * z[i - 1]) / l[i];
        }

        // Step 10
        l[N] = aArray[N] - cArray[N] * u[N - 1];
        z[N] = (dArray[N] - cArray[N] * z[N - 1]) / l[N];

        // Step 11
        const v = new Array(N + 1);
        v[N] = z[N];
        w[N] += v[N];

        // Step 12
        for (let i = N - 1; i >= 1; i--) {
            v[i] = z[i] - u[i] * v[i + 1];
            w[i] += v[i];
        }

        // Step 13
        const normV = Math.max(...v.slice(1).map(Math.abs));
        if (normV <= TOL) {
            // Success
            for (let i = 0; i <= N + 1; i++) {
                const xi = a + i * h;
                const row = document.createElement('tr');
                row.innerHTML = `
                    <td>${i}</td>
                    <td>${xi.toFixed(decimalPlaces)}</td>
                    <td>${w[i].toFixed(decimalPlaces)}</td>
                `;
                resultTable.appendChild(row);
            }
            calculationsDiv.innerHTML = `<div class="step">Success in ${k} iterations (norm v = ${normV.toExponential(2)} &lt; TOL).</div>`;
            return;
        }

        // Increment iteration
        k++;
    }

    // Step 17: If M iterations exceeded
    calculationsDiv.innerHTML = `<div class="step">Maximum number of iterations exceeded.</div>`;
}

// Helper to approximate fy = partial derivative w.r.t y
function fy(f, x, y, yPrime) {
    const h = 1e-6;
    return (f(x, y + h, yPrime) - f(x, y, yPrime)) / h;
}

// Adjust input width as user types
function adjustWidth(input) {
    input.style.width = ((input.value.length + 1) * 8) + 'px';
}
