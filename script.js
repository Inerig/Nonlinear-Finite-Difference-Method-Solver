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
            l[i] = aCoeff[i] - cCoeff[i] * u[i - 1];
            u[i] = bCoeff[i] / l[i];
            z[i] = (dCoeff[i] - cCoeff[i] * z[i - 1]) / l[i];
        }

        l[N] = aCoeff[N] - cCoeff[N] * u[N - 1];
        z[N] = (dCoeff[N] - cCoeff[N] * z[N - 1]) / l[N];

        v[N] = z[N];
        for (let i = N - 1; i >= 1; i--) {
            v[i] = z[i] - u[i] * v[i + 1];
        }

        const vNorm = Math.max(...v.slice(1).map(Math.abs));
        for (let i = 1; i <= N; i++) {
            w[i] += v[i];
        }

        calculationsDiv.innerHTML += `<div>vNorm = ${vNorm.toFixed(decimals)}</div>`;

        if (vNorm <= TOL) {
            success = true;
            break;
        }

        k++;
    }

    // Display results
    for (let i = 0; i <= N + 1; i++) {
        const row = document.createElement('tr');
        row.innerHTML = `<td>${i}</td><td>${x[i].toFixed(decimals)}</td><td>${w[i].toFixed(decimals)}</td>`;
        tbody.appendChild(row);
    }

    if (!success) {
        calculationsDiv.innerHTML += `<div style="color:red;">Maximum number of iterations exceeded!</div>`;
    } else {
        calculationsDiv.innerHTML += `<div style="color:green;">Solution found in ${k} iterations!</div>`;
    }
}
