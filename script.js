function parseFunction(funcStr) {
    try {
        // Validate the function string
        if (!/^[a-zA-Z0-9\s\+\-\*\/\(\)\.\,\[\]\^\%\&\|\!\<\>\=\~\`\@\#\$\^\\]+$/.test(funcStr)) {
            throw new Error("Invalid characters in function.");
        }
        // Create a safe function from the string input with predefined Math functions
        const mathFunctions = {
            Math: Math,
            exp: Math.exp,
            log: Math.log,
            sqrt: Math.sqrt,
            pow: Math.pow,
            sin: Math.sin,
            cos: Math.cos,
            tan: Math.tan,
            asin: Math.asin,
            acos: Math.acos,
            atan: Math.atan,
            sinh: Math.sinh,
            cosh: Math.cosh,
            tanh: Math.tanh,
            abs: Math.abs,
            floor: Math.floor,
            ceil: Math.ceil,
            round: Math.round,
            trunc: Math.trunc,
            random: Math.random
        };

        const functionString = `with(Math){return (${funcStr})}`;
        return new Function('x', 'y', 'yPrime', functionString);
    } catch (error) {
        console.error("Error parsing function:", error);
        throw new Error("Invalid function input. Please make sure your function is correctly formatted.");
    }
}

function solveAndDisplay() {
    const a = parseFloat(document.getElementById('a').value);
    const b = parseFloat(document.getElementById('b').value);
    const n = parseInt(document.getElementById('n').value);
    const alpha = parseFloat(document.getElementById('alpha').value);
    const beta = parseFloat(document.getElementById('beta').value);
    const tolerance = parseFloat(document.getElementById('tolerance').value);
    const maxIterations = parseInt(document.getElementById('maxIterations').value);
    const funcStr = document.getElementById('functionInput').value.trim();
    const decimalPlaces = parseInt(document.getElementById('decimalPlaces').value);

    if (!funcStr || isNaN(decimalPlaces) || decimalPlaces < 0 || isNaN(tolerance) || tolerance <= 0 || isNaN(maxIterations) || maxIterations <= 0) {
        alert("Please enter valid functions, non-negative integers for decimal places, positive tolerance, and a positive integer for maximum iterations.");
        return;
    }

    try {
        console.log("Function String:", funcStr);  // Debug: Check the input function string
        const f = parseFunction(funcStr);
        console.log("Parsed Function:", f.toString());  // Debug: Check the parsed function

        // Compute step size
        const h = (b - a) / (n + 1);
        // Initialize arrays to store results
        let x = new Array(n + 2);
        let w = new Array(n + 2);
        let v = new Array(n + 2);
        // Set initial conditions
        x[0] = a;
        x[n + 1] = b;
        w[0] = alpha;
        w[n + 1] = beta;

        // Set initial guesses for w_i
        for (let i = 1; i < n; i++) {
            x[i] = a + i * h;
            w[i] = alpha + (i * (beta - alpha)) / (b - a);
        }

        let k = 1;
        let converged = false;

        while (k <= maxIterations && !converged) {
            // Set coefficients for the first point
            let a1 = 2 + h * h * f(x[1], w[1], (w[2] - w[0]) / (2 * h));
            let b1 = -1 + (h / 2) * f(x[1], w[1], (w[2] - w[0]) / (2 * h));
            let d1 = -(2 * w[1] - w[2] - w[0] + h * h * f(x[1], w[1], (w[2] - w[0]) / (2 * h))) / 2;

            // Set coefficients for internal points
            for (let i = 2; i < n; i++) {
                x[i] = a + i * h;
                let t = (w[i + 1] - w[i - 1]) / (2 * h);
                let ai = 2 + h * h * f(x[i], w[i], t);
                let bi = -1 + (h / 2) * f(x[i], w[i], t);
                let ci = -1 - (h / 2) * f(x[i], w[i], t);
                let di = -(2 * w[i] - w[i + 1] - w[i - 1] + h * h * f(x[i], w[i], t)) / 2;
                // Store coefficients in arrays
                a[i] = ai;
                b[i] = bi;
                c[i] = ci;
                d[i] = di;
            }

            // Set coefficients for the last point
            x[n] = b - h;
            let t = (beta - w[n - 1]) / (2 * h);
            let aN = 2 + h * h * f(x[n], w[n], t);
            let cN = -1 - (h / 2) * f(x[n], w[n], t);
            let dN = -(2 * w[n] - w[n - 1] - beta + h * h * f(x[n], w[n], t)) / 2;

            // Solve the tridiagonal system using Thomas algorithm
            let l = new Array(n);
            let u = new Array(n);
            let z = new Array(n + 1);
            l[0] = a1;
            u[0] = b1 / a1;
            z[0] = d1 / l[0];

            for (let i = 2; i < n; i++) {
                l[i] = a[i] - c[i] * u[i - 1];
                u[i] = b[i] / l[i];
                z[i] = (d[i] - c[i] * z[i - 1]) / l[i];
            }

            l[n] = aN - cN * u[n - 1];
            z[n] = (dN - cN * z[n - 1]) / l[n];

            v[n] = z[n];

            for (let i = n - 1; i > 0; i--) {
                v[i] = z[i] - u[i] * v[i + 1];
            }

            // Update w_i
            for (let i = 0; i <= n + 1; i++) {
                w[i] = w[i] + v[i];
            }

            // Check convergence
            let maxAbsV = Math.max(...Math.abs(v));

            if (maxAbsV < tolerance) {
                converged = true;
            } else {
                k++;
            }
        }

        if (converged) {
            // Display calculations for Nonlinear Finite-Difference Method
            const nldfCalculationsContent = document.getElementById('nldfCalculationsContent');
            nldfCalculationsContent.innerHTML = '';

            for (let i = 0; i <= n + 1; i++) {
                const step = document.createElement('div');
                step.className = 'step';
                step.innerHTML = `
                    <div class="formula">Step ${i}: x = ${x[i].toFixed(decimalPlaces)}, ω = ${w[i].toFixed(decimalPlaces)}</div>
                    <div class="formula">h = ${h.toFixed(decimalPlaces)}</div>
                    <div class="formula">Coefficients:</div>
                    <div class="formula">a${i} = 2 + h²f(x${i}, w${i}, (w${i+1} - w${i-1})/(2h)) = 2 + ${h.toFixed(decimalPlaces)}² * ${f(x[i], w[i], (w[i + 1] - w[i - 1]) / (2 * h)).toFixed(decimalPlaces)} = ${a[i].toFixed(decimalPlaces)}</div>
                    <div class="formula">b${i} = -1 + (h/2)f(x${i}, w${i}, (w${i+1} - w${i-1})/(2h)) = -1 + (0.5 * ${h.toFixed(decimalPlaces)}) * ${f(x[i], w[i], (w[i + 1] - w[i - 1]) / (2 * h)).toFixed(decimalPlaces)} = ${b[i].toFixed(decimalPlaces)}</div>
                    <div class="formula">c${i} = -1 - (h/2)f(x${i}, w${i}, (w${i+1} - w${i-1})/(2h)) = -1 - (0.5 * ${h.toFixed(decimalPlaces)}) * ${f(x[i], w[i], (w[i + 1] - w[i - 1]) / (2 * h)).toFixed(decimalPlaces)} = ${c[i].toFixed(decimalPlaces)}</div>
                    <div class="formula">d${i} = -(2w${i} - w${i+1} - w${i-1} + h²f(x${i}, w${i}, (w${i+1} - w${i-1})/(2h)))/2 = -(${2 * w[i].toFixed(decimalPlaces)} - ${w[i + 1].toFixed(decimalPlaces)} - ${w[i - 1].toFixed(decimalPlaces)} + ${h.toFixed(decimalPlaces)}² * ${f(x[i], w[i], (w[i + 1] - w[i - 1]) / (2 * h)).toFixed(decimalPlaces)})/2 = ${d[i].toFixed(decimalPlaces)}</div>
                    <div class="formula">l${i} = a${i} - c${i}u${i-1} = ${a[i].toFixed(decimalPlaces)} - ${c[i].toFixed(decimalPlaces)} * ${u[i-1].toFixed(decimalPlaces)} = ${l[i].toFixed(decimalPlaces)}</div>
                    <div class="formula">u${i} = b${i}/l${i} = ${b[i].toFixed(decimalPlaces)} / ${l[i].toFixed(decimalPlaces)} = ${u[i].toFixed(decimalPlaces)}</div>
                    <div class="formula">z${i} = (d${i} - c${i}z${i-1})/l${i} = (${d[i].toFixed(decimalPlaces)} - ${c[i].toFixed(decimalPlaces)} * ${z[i-1].toFixed(decimalPlaces)}) / ${l[i].toFixed(decimalPlaces)} = ${z[i].toFixed(decimalPlaces)}</div>
                    <div class="formula">v${i} = z${i} - u${i}v${i+1} = ${z[i].toFixed(decimalPlaces)} - ${u[i].toFixed(decimalPlaces)} * ${v[i + 1].toFixed(decimalPlaces)} = ${v[i].toFixed(decimalPlaces)}</div>
                    <div class="formula">ω${i} = ω${i} + v${i} = ${w[i].toFixed(decimalPlaces)} + ${v[i].toFixed(decimalPlaces)} = ${w[i].toFixed(decimalPlaces)}</div>
                `;
                nldfCalculationsContent.appendChild(step);
            }

            // Populate the table for Nonlinear Finite-Difference Method
            const nldfTableBody = document.querySelector('#nldfResultTable tbody');
            nldfTableBody.innerHTML = '';
            for (let i = 0; i <= n + 1; i++) {
                nldfTableBody.innerHTML += `<tr><td>${i}</td><td>${x[i].toFixed(decimalPlaces)}</td><td>${w[i].toFixed(decimalPlaces)}</td></tr>`;
            }
        } else {
            alert("Maximum number of iterations exceeded.");
        }
    } catch (e) {
        alert(`Error: ${e.message}`);
    }
}

function adjustWidth(input) {
    // Calculate the width based on the number of characters
    const charCount = input.value.length;
    const minWidth = 100;
    const maxWidth = 300;
    const stepSize = (maxWidth - minWidth) / 100; // Assuming maximum character count is 100

    const newWidth = Math.min(maxWidth, minWidth + charCount * stepSize);
    input.style.width = `${newWidth}px`;
}
