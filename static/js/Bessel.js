// stricter parsing and error handling
"use strict";





//-------------------------------------------------------------------------------------------------
function CreateSymmetricMatrix(n, real = true) {
    // Create a random symmetric (real) or Hermitian (complex) matrix
    //
    // Parameters:
    //   n    : Matrix size (n×n)
    //   real : true for real symmetric, false for complex Hermitian (default: true)
    //
    // Returns:
    //   n×n symmetric (real) or Hermitian (complex) matrix
    //   Random values in range [-1, 1]
    //
    // Properties:
    //   - Symmetric: A[i,j] = A[j,i] (for real)
    //   - Hermitian: A[i,j] = conj(A[j,i]) (for complex)
    //   - Diagonal elements are always real
    //
    // Examples:
    //   CreateSymmetricMatrix(3)        → 3×3 real symmetric
    //   CreateSymmetricMatrix(4, true)  → 4×4 real symmetric
    //   CreateSymmetricMatrix(5, false) → 5×5 complex Hermitian
    //
    // Author   : Dr. Yavuz Kaya, P.Eng.
    // Modified : 16.Feb.2026

    // INPUT VALIDATION
    if (!Number.isInteger(n) || n <= 0) { throw new Error('CreateSymmetricMatrix: n must be a positive integer'); }
    
    if (typeof real !== 'boolean') { throw new Error('CreateSymmetricMatrix: real must be true or false'); }
    
    // HELPER: Random value in [-1, 1]
    function randomValue() { return 2 * Math.random() - 1; }
    
    // CREATE MATRIX
    const A = new Array(n);
    for (let i = 0; i < n; i++) {
        A[i] = new Array(n);
    }
    
    // FILL MATRIX
    if (real) {
        // REAL SYMMETRIC MATRIX
        
        // Fill upper triangle (including diagonal)
        for (let i = 0; i < n; i++) {
            for (let j = i; j < n; j++) {
                A[i][j] = randomValue();
            }
        }
        
        // Mirror to lower triangle: A[j,i] = A[i,j]
        for (let i = 0; i < n; i++) {
            for (let j = 0; j < i; j++) {
                A[i][j] = A[j][i];
            }
        }
    } else {
        // COMPLEX HERMITIAN MATRIX
        
        // Fill upper triangle (including diagonal)
        for (let i = 0; i < n; i++) {
            for (let j = i; j < n; j++) {
                if (i === j) {
                    // Diagonal: real values only
                    A[i][j] = new ComplexNum(randomValue(), 0);
                } else {
                    // Off-diagonal: complex values
                    A[i][j] = new ComplexNum(randomValue(), randomValue());
                }
            }
        }
        
        // Mirror to lower triangle: A[j,i] = conj(A[i,j])
        for (let i = 0; i < n; i++) {
            for (let j = 0; j < i; j++) {
                A[i][j] = Conj(A[j][i]);
            }
        }
    }
    
    return A;
}
//-------------------------------------------------------------------------------------------------
function randomBoolean() { return Math.random() < 0.5; }
//-------------------------------------------------------------------------------------------------
function QR_Test() {

    let maxErr, err1, err2, Ver1, Ver2, Ver3, qrRes, hesRes, isReal, isTridiagonal, isSymmHermi, ecoSize, A=[];
    let isUT, Flag, Options;
    let sum1=0, sum2=0, sum3=0, sum4=0, sum5=0, sum6=0, sum7=0, sum8=0;
    let startTime, totalTime=0, mSize = 15, tol=1e-12;
    let NumSim = 1000;

    for (let rr=0; rr<NumSim; rr++) {

        console.log(rr)
        isReal        = randomBoolean();
        isTridiagonal = randomBoolean();
        isSymmHermi   = randomBoolean();
        ecoSize       = randomBoolean();

        A         = CreateSymmetricMatrix(mSize, isReal);  
        startTime = performance.now();
        if (isReal) {
            // Real-Valued Matrix
            if (isSymmHermi) { 

                if (isTridiagonal) {
                    
                    hesRes = Hess(Copy(A),  {structure: 'symmetric', checkMatrix: false, tol:1e-10 } );

                    Options = { auto:false, isReal: isReal, isSymmHermi: isSymmHermi, isTridiagonal:isTridiagonal, ecoSize:false, tol:tol};
                    if (randomBoolean()) { qrRes  = QR(hesRes.H); } else {  qrRes  = QR(hesRes.H, Options); }
                   
                    
                    // Validation
                    Ver1   = Subtract(hesRes.H,  Multiply(qrRes.Q, qrRes.R));                       // Norm( A - QxR  )   < 1e-10
                    Ver2   = Subtract(Multiply(Transpose(qrRes.Q), qrRes.Q), Eye(qrRes.Q.length));  // Norm( Q'xQ - I )   < 1e-10
                    Ver3   = det(qrRes.Q);                                                          // Abs ( qrRes.Q  )   == 1.0
                    isUT   = isUpperTriangular(qrRes.R, 1e-10);                                     // R is upper triangular (all below-diagonal entries ≈ 0)
                    Flag   = 'QR_real_symmetric_tridiagonal';

                    qrRes.A     = hesRes.H;
                    qrRes.Flag  = Flag;
                    sum1++;
                    
                }
                else {

                    Options = { auto:false, isReal: isReal, isSymmHermi: isSymmHermi, isTridiagonal:isTridiagonal, ecoSize:false, tol:tol};
                    if (randomBoolean()) { qrRes  = QR(A); } else {  qrRes  = QR(A, Options); }

                    // Validation
                    Ver1   = Subtract(A,  Multiply(qrRes.Q, qrRes.R));                                // Norm( A - QxR  )   < 1e-10
                    Ver2   = Subtract(Multiply(Transpose(qrRes.Q), qrRes.Q), Eye(qrRes.Q.length));    // Norm( Q'xQ - I )   < 1e-10
                    Ver3   = det(qrRes.Q);                                                            // Abs ( qrRes.Q  )   == 1.0
                    isUT   = isUpperTriangular(qrRes.R, 1e-10);                                       // R is upper triangular (all below-diagonal entries ≈ 0)
                    Flag   = 'QR_real_symmetric';
                    
                    qrRes.A     = A;
                    qrRes.Flag  = Flag;
                    sum2++;
                }
            
            } 
            else {
                
                A[2][0] = 10.54;  // non-symmetric 

                if (ecoSize) {

                    if (randomBoolean()) { 
                        // tall rectangular matrix  (add two rows)
                        A.push(Rand(mSize)); A.push(Rand(mSize));
                        Flag   = 'QR_real_eco_tall';
                    } else { 
                        // wide rectangular matrix
                        for (let i = 0; i < A.length; i++) { A[i].push(Math.random(), Math.random());}
                        Flag   = 'QR_real_full';
                    } 
                    
                    Options = { auto:false, isReal: isReal, isSymmHermi: isSymmHermi, isTridiagonal:false, ecoSize:ecoSize, tol:tol};
                    if (randomBoolean()) { qrRes  = QR(A); } else {  qrRes  = QR(A, Options); }
                   
                    // Validation
                    Ver1   = Subtract(A,  Multiply(qrRes.Q, qrRes.R));                               // Norm( A - QxR  )   < 1e-10
                    Ver2   = Subtract(Multiply(Transpose(qrRes.Q), qrRes.Q), Eye(qrRes.Q.length));   // Norm( Q'xQ - I )   < 1e-10
                    if (qrRes.Q.length != qrRes.Q[0].length) {Ver3 = 1; } else { Ver3   = det(qrRes.Q); } // Abs ( qrRes.Q  )   == 1.0
                    isUT   = isUpperTriangular(qrRes.R, 1e-10);                                      // R is upper triangular (all below-diagonal entries ≈ 0)

                    qrRes.A     = A;
                    qrRes.Flag  = Flag;
                    sum3++;
                }
                else {

                    Options = { auto:false, isReal: isReal, isSymmHermi: isSymmHermi, isTridiagonal:false, ecoSize:ecoSize, tol:tol};
                    if (randomBoolean()) { qrRes  = QR(A); } else {  qrRes  = QR(A, Options); }

                    // Validation
                    Ver1   = Subtract(A,  Multiply(qrRes.Q, qrRes.R));                               // Norm( A - QxR  )   < 1e-10
                    Ver2   = Subtract(Multiply(Transpose(qrRes.Q), qrRes.Q), Eye(qrRes.Q.length));   // Norm( Q'xQ - I )   < 1e-10
                    Ver3   = det(qrRes.Q);                                                           // Abs ( qrRes.Q  )   == 1.0
                    isUT   = isUpperTriangular(qrRes.R, 1e-10);                                      // R is upper triangular (all below-diagonal entries ≈ 0)
                    Flag   = 'QR_real_full';
                    
                    qrRes.A     = A;
                    qrRes.Flag  = Flag;
                    sum4++;
                }

            }
        } else {

            // Complex-Valued Matrix
            if (isSymmHermi) {

                if (isTridiagonal) {

                    hesRes = Hess(Copy(A),  {structure: 'hermitian', checkMatrix: false, tol:1e-10 } );
                    Options = { auto:false, isReal: isReal, isSymmHermi: isSymmHermi, isTridiagonal:isTridiagonal, ecoSize:false, tol:tol};
                    if (randomBoolean()) { qrRes  = QR(hesRes.H); } else {  qrRes  = QR(hesRes.H, Options); }

                    // Validation
                    Ver1   = Subtract(hesRes.H,  Multiply(qrRes.Q, qrRes.R));                       // Norm( A - QxR  )   < 1e-10
                    Ver2   = Subtract(Multiply(Transpose(qrRes.Q), qrRes.Q), Eye(qrRes.Q.length));  // Norm( Q'xQ - I )   < 1e-10
                    Ver3   = det(qrRes.Q);                                                          // Abs ( qrRes.Q  )   == 1.0
                    isUT   = isUpperTriangular(qrRes.R, 1e-10);                                     // R is upper triangular (all below-diagonal entries ≈ 0)
                    Flag   = 'QR_complex_hermitian_tridiagonal';
                    sum5++;
                }
                else {
                    Options = { auto:false, isReal: isReal, isSymmHermi: isSymmHermi, isTridiagonal:isTridiagonal, ecoSize:false, tol:tol};
                    if (randomBoolean()) { qrRes  = QR(A); } else {  qrRes  = QR(A, Options); }
                    
                    // Validation
                    Ver1   = Subtract(A,  Multiply(qrRes.Q, qrRes.R));                              // Norm( A - QxR  )   < 1e-10
                    Ver2   = Subtract(Multiply(Transpose(qrRes.Q), qrRes.Q), Eye(qrRes.Q.length));  // Norm( Q'xQ - I )   < 1e-10
                    Ver3   = det(qrRes.Q);                                                          // Abs ( qrRes.Q  )   == 1.0
                    isUT   = isUpperTriangular(qrRes.R, 1e-10);                                     // R is upper triangular (all below-diagonal entries ≈ 0)
                    Flag   = 'QR_complex_hermitian';

                    qrRes.A     = A;
                    qrRes.Flag  = Flag;
                    sum6++;
                }
                
            } 
            else {
                
                A[2][0] = new ComplexNum(10.54, 8.98);  // non-symmetric 

                if (ecoSize) {

                    if (randomBoolean()) { 
                        // tall rectangular matrix  (add two rows)
                        A.push(Rand_Complex(mSize)); A.push(Rand_Complex(mSize));
                        Flag   = 'QR_complex_eco_tall';
                    } else { 
                        // wide rectangular matrix
                        for (let i = 0; i < A.length; i++) { A[i].push(Rand_Complex(1)[0], Rand_Complex(1)[0]); }  
                        Flag   = 'QR_complex_full';
                    }
                    
                    Options = { auto:false, isReal: isReal, isSymmHermi: isSymmHermi, isTridiagonal:false, ecoSize:ecoSize, tol:tol};
                    if (randomBoolean()) { qrRes  = QR(A); } else {  qrRes  = QR(A, Options); }

                    // Validation
                    Ver1   = Subtract(A,  Multiply(qrRes.Q, qrRes.R));                               // Norm( A - QxR  )   < 1e-10
                    Ver2   = Subtract(Multiply(Transpose(qrRes.Q), qrRes.Q), Eye(qrRes.Q.length));   // Norm( Q'xQ - I )   < 1e-10
                    if (qrRes.Q.length != qrRes.Q[0].length) {Ver3 = 1; } else { Ver3   = det(qrRes.Q); } // Abs ( qrRes.Q  )   == 1.0
                    isUT   = isUpperTriangular(qrRes.R, 1e-10);                                      // R is upper triangular (all below-diagonal entries ≈ 0)

                    qrRes.A     = A;
                    qrRes.Flag  = Flag;
                    sum7++;

                }
                else {

                    Options = { auto:false, isReal: isReal, isSymmHermi: isSymmHermi, isTridiagonal:false, ecoSize:ecoSize, tol:tol};
                    if (randomBoolean()) { qrRes  = QR(A); } else {  qrRes  = QR(A, Options); }

                    // Validation
                    Ver1   = Subtract(A,  Multiply(qrRes.Q, qrRes.R));                               // Norm( A - QxR  )   < 1e-10
                    Ver2   = Subtract(Multiply(Transpose(qrRes.Q), qrRes.Q), Eye(qrRes.Q.length));   // Norm( Q'xQ - I )   < 1e-10
                    Ver3   = det(qrRes.Q);                                                           // Abs ( qrRes.Q  )   == 1.0
                    isUT   = isUpperTriangular(qrRes.R, 1e-10);                                      // R is upper triangular (all below-diagonal entries ≈ 0)
                    Flag   = 'QR_complex_full';

                    qrRes.A     = A;
                    qrRes.Flag  = Flag;
                    sum8++;
                }
                
            } 
        }
        totalTime += (performance.now() - startTime);

        
        // QR results
        console.log(qrRes)
        console.log(Flag)
        console.log('maxErr : ' + maxErr );
        
        // Error calculation
        err1   = Norm(Ver1, 'fro');
        err2   = Norm(Ver2, 'fro');
        maxErr = Math.max(err1, err2);

        // Check error     
        if ((maxErr > 1e-10 || !isUT || (Abs(Ver3)-1 > 1e-10) )) {

            //Print(A, !isReal)
            console.log(Flag)
            console.log('maxErr : ' + maxErr );

            break;
        }
        console.log('-------------------------------------------------------------------')
    
    }

    console.log('-------------------------------------------------------------------')
    console.log('---                    SUMMARY                                  ---')
    console.log('-------------------------------------------------------------------')
    console.log('Real-----Symmetric (Square)----------Tridiagonal--(sum1)  : ' + sum1);
    console.log('Real-----Symmetric (Square)-----------------------(sum2)  : ' + sum2);
    console.log('Real-----NonSymmetric (Rectangular)--EconomySize--(sum3)  : ' + sum3);
    console.log('Real-----NonSymmetric (Rectangular)---------------(sum4)  : ' + sum4);

    console.log('Complex--Hermitian (Square)----------Tridiagonal--(sum5)  : ' + sum5);
    console.log('Complex--Hermitian (Square)-----------------------(sum6)  : ' + sum6);
    console.log('Complex--NonHermitian (Rectangular)--EconomySize--(sum7)  : ' + sum7);
    console.log('Complex--NonHermitian (Rectangular)---------------(sum8)  : ' + sum8);

    console.log('---------------------------------------------------Total  : ' + (sum1+sum2+sum3+sum4+sum5+sum6+sum7+sum8).toString())
    console.log('Average Time (ms) : ' + (totalTime / NumSim).toPrecision(2))


    function isUpperTriangular(matrix, tol) {
        for (let i = 1; i < matrix.length; i++) {
            for (let j = 0; j < i; j++) {
                const val = (matrix[i][j] instanceof ComplexNum) ? Math.sqrt(matrix[i][j].Re**2 + matrix[i][j].Im**2) : Math.abs(matrix[i][j]);
                if (val > tol) return false;
            }
        }
        return true;
    }

}
//-------------------------------------------------------------------------------------------------
function Hess_Test() {
    let A=[], isReal, isSymmHermi, Hess_res, flag=false;
    let mSize = 37,   tol = 1e-10;
    let NumSim = 100;

    for (let rr=0; rr<NumSim; rr++) {

        isReal       = randomBoolean();
        isSymmHermi  = randomBoolean();
        A            = CreateSymmetricMatrix(mSize, isReal); 

        if (isReal) {
            // Real-Valued Matrix
            if (isSymmHermi) {
                // Real-Valued-Symmetric
                Hess_res = Hess(A, {auto: false,   isReal:isReal,   isSymmHermi:isSymmHermi,   tol:1e-14});
            }
            else {
                // Real-Valued-NonSymmetric
                A[2][0] = 10.5445;  // non-symmetric
                A[3][0] = 897.533;
                Hess_res = Hess(A, {auto: false,   isReal:isReal,   isSymmHermi:isSymmHermi,   tol:1e-14});
            }
        } else {
            // Complex-Valued Matrix
            if (isSymmHermi) {
                // Complex-Valued-Hermetian
                Hess_res = Hess(A, {auto: false,   isReal:isReal,   isSymmHermi:isSymmHermi,   tol:1e-14});
            }
            else {
                // Complex-Valued-NonHermetian
                A[2][0] = new ComplexNum(10.54, 8.98);  // non-symmetric 
                A[3][0] = new ComplexNum(879, 898);
                Hess_res = Hess(A, {auto: false,   isReal:isReal,   isSymmHermi:isSymmHermi,   tol:1e-14});
            }
        }

        // Verify the results
        let t1 = Multiply(Transpose(Hess_res.Q), Hess_res.Q);                       // I = Q^T * Q
        let t2 = Multiply(Multiply(Hess_res.Q, Hess_res.H), Transpose(Hess_res.Q))  // A = Q * H * Q^T
        let t3 = Norm(Subtract(t2, A),           'fro');                            //  ||Q*H*Q^H - A||_F
        let t4 = Norm(Subtract(t1, Eye(mSize)),  'fro');                            // ||Q^H*Q - I||_F

        if (t3 > tol  ) {  flag = true;  }
        if (t4 > tol  ) {  flag = true;  }

        console.log('Iteration Number   : ' + rr)
        if (flag) {
            console.log('Iteration Number      : ' + rr)
            console.log('IsReal                : ' + isReal)
            console.log('isSymmHermi           : ' + isSymmHermi)
            console.log('I = Q^T * Q      err  : ' + t4)
            console.log('A = Q * H * Q^T  err  : ' + t3)
            return
        }
    }
    return "ASD"
}
//-------------------------------------------------------------------------------------------------
function det(A) {
    const n = A.length;

    // Validate square matrix
    if (A.some(row => row.length !== n)) {
        throw new Error("Matrix must be square");
    }

    // Normalize entries: wrap plain numbers as ComplexNum
    const toC = v => (v instanceof ComplexNum) ? v : new ComplexNum(v, 0);

    // Work on a deep copy so we don't mutate the original
    let M = A.map(row => row.map(toC));

    // LU decomposition with partial pivoting (Gaussian elimination)
    // Track sign flips from row swaps
    let sign = new ComplexNum(1, 0);

    for (let col = 0; col < n; col++) {

        // --- Partial pivoting: find row with largest |pivot| at or below current row ---
        let maxAbs = -1, pivotRow = -1;
        for (let row = col; row < n; row++) {
            const a = M[row][col].Abs();
            if (a > maxAbs) { maxAbs = a; pivotRow = row; }
        }

        // Singular matrix
        if (maxAbs === 0) return new ComplexNum(0, 0);

        // Swap rows if needed
        if (pivotRow !== col) {
            [M[col], M[pivotRow]] = [M[pivotRow], M[col]];
            sign = sign.Multiply(-1);   // each swap flips the sign of det
        }

        const pivot = M[col][col];

        // Eliminate entries below the pivot
        for (let row = col + 1; row < n; row++) {
            const factor = M[row][col].Divide(pivot);
            for (let k = col; k < n; k++) {
                M[row][k] = M[row][k].Subtract(factor.Multiply(M[col][k]));
            }
        }
    }

    // det = sign * product of diagonal entries (upper triangular after elimination)
    let result = sign;
    for (let i = 0; i < n; i++) {
        result = result.Multiply(M[i][i]);
    }

    return result; // ComplexNum
}
//-------------------------------------------------------------------------------------------------
function Eig_Test() {

    let A=[], B=[], isReal, isSymmHermi, isPosDef, isTridiagonal, label;
    let r1=0, r2=0, r3=0, r4=0, r5=0, r6=0, r7=0, r8=0;
    let Hess_res, flag=false;
    let tol = 1e-14;
    let NumSim = 10000;
    let mMax = 5, mMin=10;

    // Counters
    let nPass = 0, nFail = 0;
    const failures = [];

    for (let rr=0; rr<NumSim; rr++) {

        // Randomly select the size of matrix
        let mSize = Math.floor(Math.random() * (mMax - mMin + 1)) + mMin;

        // Randomly choose real or complex and Build a random base matrix-B
        isReal         = randomBoolean();                                                           isReal=true;
        B              = isReal ? Create_Real_Matrix(mSize) : Create_Complex_Matrix(mSize);

        // Randomly choose matrix class; construction follows the flag
        isSymmHermi    = randomBoolean();                                                           isSymmHermi=false;
        isPosDef       = isSymmHermi && randomBoolean();    // PositiveDefinite requires isSymmHermi to be true
        isTridiagonal  = randomBoolean();
        
        // ---- Construct A to match the chosen class ----------------------
        if (isReal) {
            // Real-valued Matrix
            if (isSymmHermi) {
                A = MakeSymmetric(B);
                if(isPosDef) { A = MakePosDef(B); }

                if (isTridiagonal) {
                    Hess_res = Hess(A, {auto: false, isReal: isReal, isSymmHermi: isSymmHermi, tol: tol});
                    A = Hess_res.H;
                    r1++;
                    label = 'Real | Symmetric | Tridiagonal';
                } else {
                    if (isPosDef) {
                        r2++;
                        label = 'Real | Symmetric | PosDef';
                    }
                    else {
                        r3++;
                        label = 'Real | Symmetric';
                    }
                }
            } else {
                A = B;
                r4++;
                label = 'Real | General';
            }
        } else {
            if (isSymmHermi) {
                A = MakeSymmetric(B);
                if(isPosDef) { A = MakePosDef(B); }

                if (isTridiagonal) {
                    Hess_res = Hess(A, { auto: false, isReal: isReal, isSymmHermi: isSymmHermi, tol: tol });
                    A = Hess_res.H;
                    r5++;
                    label = 'Complex | Hermitian | Tridiagonal';
                    
                } else {
                    if (isPosDef) {
                        r6++;
                        label = 'Complex | Hermitian | PosDef';
                    }
                    else {
                        r7++;
                        label = 'Complex | Hermitian';
                    }
                }
            } else {
                A = B;
                r8++;
                label = 'Complex | General';
            }
        }

        // ---- Call Eig --------------------------------------------------
        let res = Eig(A, { auto: false, isReal: isReal, isSymmHermi: isSymmHermi, isPosDef: isPosDef, isTridiagonal: isTridiagonal, tol: tol });

        // ---- Skip unimplemented helpers (stubs return undefined) -------
        if (res == null || res.eigenvalues == null || res.eigenvectors == null) { continue; }
        
        // ---- Verification ----------------------------------------------
        let Ver_result = Verify_Eig(A, res.eigenvalues, res.eigenvectors, tol * mSize)

        console.log(rr)
        console.log('IsReal               : ' + isReal)
        console.log('isSymmHermi          : ' + isSymmHermi)
        console.log('isPosDef             : ' + isPosDef)
        console.log('isTridiagonal        : ' + isTridiagonal)
        console.log('Ver_result           : ' + Ver_result.pass)

        if (Ver_result.pass == false) { 
            if (Ver_result.pass == false) { return Ver_result; }
        }

        console.log('--------------------------------------------------')
    }


    return "ASD";

    // helper functions
    function Create_Real_Matrix(m) {
        return Array.from({ length: m }, () => Array.from({ length: m }, () => (Math.random() * 2 - 1)));
    }

    function Create_Complex_Matrix(m) {
        return Array.from({ length: m }, () => Array.from({ length: m }, () => new ComplexNum(Math.random() * 2 - 1, Math.random() * 2 - 1)));
    }

    function MakeSymmetric(B) {
        return Multiply(Add(Transpose(B), B), 0.5);
    }

    function MakePosDef(B) {
        return Multiply(Transpose(B), B);
    }

    
}
function Verify_Eig(A, eigenvalues, eigenvectors, tol) {
    let D = eigenvalues;
    let V = Transpose(eigenvectors, {conjugate: false});
    let M1 = Multiply(A, V);
    let M2 = Multiply(V, Diag(D));
    let Diff = Subtract(M1, M2);
    let nEig = eigenvalues.length;
    let residual = new Array(nEig);
    let pass     = new Array(nEig);
    let allPass  = true;
    
    // Compute ||A||_F once for size-independent relative residual
    const normA = Norm(A, 'fro');
    
    for (let i = 0; i < nEig; i++) {
        let col = Diff.map(row => row[i]);
        residual[i] = Norm(col, 'fro') / normA;  // relative to matrix norm
        pass[i]     = residual[i] < tol;
        if (!pass[i]) { allPass = false; }
    }
    return {
        matrix       : A,
        pass         : allPass,
        eigenvalues  : D,
        residual     : residual,
        perEigenPass : pass
    };
}
//-------------------------------------------------------------------------------------------------
function Eig(A, Option) {

    // Computes the eigenvalues and eigenvectors of a square matrix A.
    // Automatically detects matrix properties (real/complex, symmetric/Hermitian,
    // positive definite, tridiagonal) and dispatches to the appropriate helper.
    //
    // Parameters:
    //   A      : 2D array — square matrix (real or complex-valued)
    //   Option : options object (optional)
    //     {
    //       auto          : boolean  →  true  = auto-detect all matrix properties [default: true]
    //                                   false = use caller-supplied flags below (no re-validation)
    //       isReal        : boolean  →  true  = all elements are real-valued
    //                                   false = matrix contains ComplexNum elements
    //       isSymmHermi   : boolean  →  true  = symmetric (real) or Hermitian (complex)
    //       isPosDef      : boolean  →  true  = symmetric/Hermitian and positive definite
    //       isTridiagonal : boolean  →  true  = matrix is tridiagonal
    //       nEig          : integer  →  number of eigenpairs to compute, selected by smallest
    //                                   magnitude; omit or set >= n to compute all [default: n]
    //       tol           : number   →  numerical zero threshold  [default: 1e-14]
    //     }
    //
    // Returns:  { Val, Vec }
    //   Val : 1D array of eigenvalues (real numbers or ComplexNum objects), length nEig
    //   Vec : 2D array whose columns are the corresponding eigenvectors, width nEig
    //
    // nEig behaviour:
    //   - Omitted or >= n  → all n eigenpairs are computed and returned
    //   - 1 <= nEig < n    → each helper computes only nEig eigenpairs internally,
    //                         returning the nEig smallest-magnitude pairs sorted in
    //                         ascending order of |eigenvalue|
    //   - nEig < 1         → throws RangeError
    //
    // Dispatch hierarchy (evaluated in order, first match wins):
    //   n === 1                        → EigScalar                      (trivial)
    //   n === 2                        → Eig2x2                         (closed-form)
    //   isReal
    //     isSymmHermi && isTridiagonal → Eig_Real_Symm_Tridiag          (direct QR on tridiagonal)
    //     isSymmHermi && isPosDef      → Eig_Real_Symm_PosDef           (Cholesky-based)
    //     isSymmHermi                  → Eig_Real_Symm                  (Hess → tridiagonal → QR)
    //                                  → Eig_Real_General               (Hess → Francis double-shift QR)
    //   isComplex
    //     isSymmHermi && isTridiagonal → Eig_Complex_Hermitian_Tridiag  (direct QR on tridiagonal)
    //     isSymmHermi && isPosDef      → Eig_Complex_Hermitian_PosDef   (Cholesky-based)
    //     isSymmHermi                  → Eig_Complex_Hermitian          (Hess → tridiagonal → QR)
    //                                  → Eig_Complex_General            (Hess → single-shift complex QR)
    //
    // Examples:
    //   Eig([[2, 1], [1, 2]])                               → { Val: [1, 3], Vec: [[-0.707, 0.707], [0.707, 0.707]] }
    //   Eig([[4, 1, 0], [1, 3, 1], [0, 1, 2]], { nEig: 2 }) → { Val: [λ₁, λ₂], Vec: [[...], [...]] }  — 2 smallest-magnitude eigenpairs
    //   Eig(A, { auto: false, isReal: true, isSymmHermi: true, isPosDef: false, isTridiagonal: false, nEig: 3 })
    //       → routes directly to Eig_Real_Symm, computes and returns 3 smallest-magnitude eigenpairs
    //
    // Author   : Dr. Yavuz Kaya, P.Eng.
    // Modified : 15.Jun.2026

    // -------------------------------------------------------------------------
    // Step 1 — Apply defaults if Option is omitted
    // -------------------------------------------------------------------------
    if (Option == null) { return Eig(A, { auto: true, tol: 1e-14 }); }

    // -------------------------------------------------------------------------
    // Step 2 — Validate Option type
    // -------------------------------------------------------------------------
    if (typeof Option !== 'object' || Array.isArray(Option)) { throw new TypeError('Eig: Option must be a plain object.'); }

    // -------------------------------------------------------------------------
    // Step 3 — Validate input matrix
    // -------------------------------------------------------------------------
    if (!Array.isArray(A) || !Array.isArray(A[0])) { throw new Error('Eig: A must be a 2D array.');      }
    if (A.length !== A[0].length)                  { throw new Error('Eig: A must be a square matrix.'); }

    // -------------------------------------------------------------------------
    // Step 4 — Resolve options
    // -------------------------------------------------------------------------
    const tol  = (Option.tol  != null) ? Option.tol  : 1e-14;
    const auto = (Option.auto != null) ? Option.auto : true;

    const n = A.length;

    // Resolve and validate nEig; clamp to n if caller passes a value >= n
    const nEig = (Option.nEig != null) ? Option.nEig : n;
    if (!Number.isInteger(nEig) || nEig < 1) { throw new RangeError(`Eig: nEig must be a positive integer (got ${nEig}).`); }
    const nEigResolved = Math.min(nEig, n);   // passed to every helper; helpers never see a value > n

    let isReal, isSymmHermi, isPosDef, isTridiagonal;

    if (auto) {
        // Auto-detect all matrix properties
        isReal        = IsReal(A);
        isSymmHermi   = isReal ? IsSymmetric(A) : IsHermitian(A);
        isPosDef      = isSymmHermi && IsPosDef(A);
        isTridiagonal = IsTridiagonal(A, tol, !isReal);
    } else {
        // Trust caller-supplied flags — no re-validation
        isReal        = (Option.isReal        != null) ? Option.isReal        : IsReal(A);
        isSymmHermi   = (Option.isSymmHermi   != null) ? Option.isSymmHermi   : (isReal ? IsSymmetric(A) : IsHermitian(A));
        isPosDef      = (Option.isPosDef      != null) ? Option.isPosDef      : false;
        isTridiagonal = (Option.isTridiagonal != null) ? Option.isTridiagonal : false;
    }

    // -------------------------------------------------------------------------
    // Step 5 — Dispatch  (nEigResolved is forwarded to every helper)
    // -------------------------------------------------------------------------

    // Scalar (1×1)
    if (n === 1) { return EigScalar(A); }

    // 2×2 — closed-form solution (real or complex)
    if (n === 2) { return Eig2x2(A, isReal, nEigResolved); }

    if (isReal) {

        
        if      (isSymmHermi && isTridiagonal) { return Eig_Real_Symm_Tridiag(A, nEigResolved);  }  // Real-Symmetric-Tridiagonal
        else if (isSymmHermi && isPosDef)      { return Eig_Real_Symm_PosDef(A, nEigResolved);   }  // Real-Symmetric_PositiveDefinite
        else if (isSymmHermi)                  { return Eig_Real_Symm_PosDef(A, nEigResolved);   }  // Real symmetric ==> same as Real-Symmetric_PositiveDefinite
        else                                   { console.log("Eig_Real_General"); 
                                                 return Eig_Real_General(A, nEigResolved, tol);       

                                                
                                                }  // Real General -- Algorithm already includes Hessenberg reduction 

    } else {

        // Complex Hermitian tridiagonal — direct QR iteration, no Hess reduction needed
        if      (isSymmHermi && isTridiagonal) { return Eig_Complex_Hermitian_Tridiag(A, nEigResolved, tol); }

        // Complex Hermitian positive definite — Cholesky-based approach
        else if (isSymmHermi && isPosDef)      { return Eig_Complex_Hermitian_PosDef(A, nEigResolved, tol);  }

        // Complex Hermitian (general) — Hess reduces to tridiagonal, then QR iteration
        else if (isSymmHermi)                  { return Eig_Complex_Hermitian(A, nEigResolved, tol);         }

        // Complex general — Hess to upper Hessenberg, then single-shift complex QR
        else                                   { return Eig_Complex_General(A, nEigResolved, tol);           }

    }

    // Each helper receives nEig and is responsible for:
    //   1. Computing only the nEig smallest-magnitude eigenpairs internally
    //   2. Returning { Val, Vec } with exactly nEig entries, sorted ascending by |Val[i]|
    
    function Eig_Real_Symm_Tridiag(A, nEig) {

        const n = A.length;
        const eps = Number.EPSILON;    // ≈ 2.22e-16 

        /* ================================================================
        1.  Extract diagonal (d) and sub-diagonal (e) vectors
        ================================================================ */
        const d = new Array(n);          // d[i] = A[i][i]
        const e = new Array(n);          // e[i] = A[i][i+1],  e[n-1] = 0

        for (let i = 0; i < n; i++) {
            d[i] = A[i][i];
            e[i] = i < n - 1 ? A[i][i + 1] : 0.0;
        }

        /* ================================================================
        2.  Eigenvector accumulator  Z  ←  Iₙ
            Column j of Z will converge to the eigenvector for d[j].
        ================================================================ */
        const Z = new Array(n);
        for (let i = 0; i < n; i++) {
            Z[i] = new Array(n).fill(0.0);
            Z[i][i] = 1.0;
        }

        /* ================================================================
        3.  Implicit QL iterations with Wilkinson shift
        ================================================================ */
        const maxIter = 30 * n;          // generous per-eigenvalue budget

        for (let l = 0; l < n; l++) {    // deflate eigenvalue at position l
            let iter = 0;

            // --- outer convergence loop for eigenvalue l ---
            while (true) {

                // 3a. find smallest m ≥ l where |e[m]| is negligible
                let m = l;
                while (m < n - 1) {
                    const offDiagTest = Math.abs(d[m]) + Math.abs(d[m + 1]);
                    if (Math.abs(e[m]) <= eps * offDiagTest) break;
                    m++;
                }

                // eigenvalue d[l] has converged
                if (m === l) break;

                if (++iter > maxIter) {
                    throw new Error(
                        "Eig_Real_Symm_Tridiag: no convergence after "
                        + maxIter + " iterations (index " + l + ")"
                    );
                }

                // 3b. Wilkinson shift  (closest eigenvalue of trailing 2×2 block)
                let g = (d[l + 1] - d[l]) / (2.0 * e[l]);
                let r = Math.hypot(g, 1.0);               // √(g² + 1)
                // shift σ  =  d[m] − d[l] + e[l] / (g + sign(g)·r)
                g = d[m] - d[l] + e[l] / (g + (g >= 0.0 ? r : -r));

                let s = 1.0;
                let c = 1.0;
                let p = 0.0;
                let underflow = false;

                // 3c. QL rotation sweep  i = m-1, m-2, … , l
                for (let i = m - 1; i >= l; i--) {
                    const f = s * e[i];
                    const b = c * e[i];
                    r = Math.hypot(f, g);
                    e[i + 1] = r;

                    // Guard: if r collapses to zero the rotation is degenerate
                    if (r === 0.0) {
                        d[i + 1] -= p;
                        e[m] = 0.0;
                        underflow = true;
                        break;                              // restart this l
                    }

                    s = f / r;
                    c = g / r;
                    g = d[i + 1] - p;
                    r = (d[i] - g) * s + 2.0 * c * b;
                    p = s * r;
                    d[i + 1] = g + p;
                    g = c * r - b;

                    // ── accumulate Givens rotation into Z ──
                    for (let k = 0; k < n; k++) {
                        const t      = Z[k][i + 1];
                        Z[k][i + 1]  = s * Z[k][i] + c * t;
                        Z[k][i]      = c * Z[k][i] - s * t;
                    }
                }

                // update diagonal & sub-diagonal (skip on underflow restart)
                if (!underflow) {
                    d[l] -= p;
                    e[l]  = g;
                    e[m]  = 0.0;
                }
            }
        }

        /* ================================================================
        4.  Sort eigenvalues in ascending order
        ================================================================ */
        const idx = Array.from({ length: n }, (_, i) => i);
        idx.sort((a, b) => d[a] - d[b]);

        /* ================================================================
        5.  Pack the first nEig eigen-pairs
        ================================================================ */
        const eigenvalues  = new Array(nEig);
        const eigenvectors = new Array(nEig);

        for (let j = 0; j < nEig; j++) {
            const col = idx[j];
            eigenvalues[j] = d[col];

            const vec = new Array(n);
            for (let k = 0; k < n; k++) vec[k] = Z[k][col];
            eigenvectors[j] = vec;
        }

        return { eigenvalues, eigenvectors };
    }
    
    function Eig_Real_Symm_PosDef(A, nEig) {

        const n   = A.length;
        const eps = Number.EPSILON;                 // ≈ 2.22e-16

        /* ================================================================
        1.  Deep-copy A → H  (H is overwritten during reduction)
        ================================================================ */
        const H = [];
        for (let i = 0; i < n; i++) {
            H[i] = new Array(n);
            for (let j = 0; j < n; j++) H[i][j] = A[i][j];
        }

        /* ================================================================
        2.  Q ← Iₙ   (accumulates Householder + Givens transforms)
            Final columns of Q are eigenvectors of A.
        ================================================================ */
        const Q = [];
        for (let i = 0; i < n; i++) {
            Q[i] = new Array(n).fill(0.0);
            Q[i][i] = 1.0;
        }

        /* ================================================================
        3.  Householder tridiagonalization   A = Q·T·Qᵀ
            Each step k zeroes H[k+2:n, k] via  H ← Pₖ·H·Pₖ
            with reflector  Pₖ = I − 2vvᵀ  acting on rows/cols k+1…n−1.
        ================================================================ */
        for (let k = 0; k < n - 2; k++) {

            const m = n - k - 1;                   // reflector dimension

            // ── extract x = H[k+1:n, k] ──
            const x = new Array(m);
            for (let i = 0; i < m; i++) x[i] = H[k + 1 + i][k];

            // ── ‖x‖₂ ──
            let xNorm = 0.0;
            for (let i = 0; i < m; i++) xNorm += x[i] * x[i];
            xNorm = Math.sqrt(xNorm);

            if (xNorm < eps) continue;              // column already zero

            // ── Householder vector v = x − α·e₁ ──
            //    sign chosen to avoid catastrophic cancellation
            const alpha = x[0] >= 0.0 ? -xNorm : xNorm;
            const v = x.slice();
            v[0] -= alpha;

            // ── normalise v ──
            let vNorm = 0.0;
            for (let i = 0; i < m; i++) vNorm += v[i] * v[i];
            vNorm = Math.sqrt(vNorm);
            if (vNorm < eps) continue;
            for (let i = 0; i < m; i++) v[i] /= vNorm;

            // ── H ← P·H   (left multiply: rows k+1…n−1) ──
            for (let j = 0; j < n; j++) {
                let dot = 0.0;
                for (let i = 0; i < m; i++) dot += v[i] * H[k + 1 + i][j];
                dot *= 2.0;
                for (let i = 0; i < m; i++) H[k + 1 + i][j] -= dot * v[i];
            }

            // ── H ← H·P   (right multiply: cols k+1…n−1) ──
            for (let i = 0; i < n; i++) {
                let dot = 0.0;
                for (let j = 0; j < m; j++) dot += H[i][k + 1 + j] * v[j];
                dot *= 2.0;
                for (let j = 0; j < m; j++) H[i][k + 1 + j] -= dot * v[j];
            }

            // ── Q ← Q·P   (accumulate transformation) ──
            for (let i = 0; i < n; i++) {
                let dot = 0.0;
                for (let j = 0; j < m; j++) dot += Q[i][k + 1 + j] * v[j];
                dot *= 2.0;
                for (let j = 0; j < m; j++) Q[i][k + 1 + j] -= dot * v[j];
            }
        }

        /* ================================================================
        4.  Extract diagonal (d) and sub-diagonal (e) from tridiagonal H
        ================================================================ */
        const d = new Array(n);
        const e = new Array(n);
        for (let i = 0; i < n; i++) {
            d[i] = H[i][i];
            e[i] = i < n - 1 ? H[i][i + 1] : 0.0;
        }

        /* ================================================================
        5.  Implicit QL iterations with Wilkinson shift
            Givens rotations are accumulated directly into Q, so
            at convergence Q holds the eigenvectors of A (= Qₕ · Z).
        ================================================================ */
        const maxIter = 30 * n;

        for (let l = 0; l < n; l++) {
            let iter = 0;

            while (true) {

                // 5a. find smallest m ≥ l where |e[m]| is negligible
                let m = l;
                while (m < n - 1) {
                    const offDiag = Math.abs(d[m]) + Math.abs(d[m + 1]);
                    if (Math.abs(e[m]) <= eps * offDiag) break;
                    m++;
                }

                if (m === l) break;                  // eigenvalue d[l] converged

                if (++iter > maxIter) {
                    throw new Error(
                        'Eig_Real_Symm_PosDef: no convergence after '
                        + maxIter + ' iterations (index ' + l + ')'
                    );
                }

                // 5b. Wilkinson shift
                let g = (d[l + 1] - d[l]) / (2.0 * e[l]);
                let r = Math.hypot(g, 1.0);
                g = d[m] - d[l] + e[l] / (g + (g >= 0.0 ? r : -r));

                let s = 1.0;
                let c = 1.0;
                let p = 0.0;
                let underflow = false;

                // 5c. QL rotation sweep  i = m−1, m−2, … , l
                for (let i = m - 1; i >= l; i--) {
                    const f = s * e[i];
                    const b = c * e[i];
                    r = Math.hypot(f, g);
                    e[i + 1] = r;

                    if (r === 0.0) {                 // degenerate rotation
                        d[i + 1] -= p;
                        e[m] = 0.0;
                        underflow = true;
                        break;
                    }

                    s = f / r;
                    c = g / r;
                    g = d[i + 1] - p;
                    r = (d[i] - g) * s + 2.0 * c * b;
                    p = s * r;
                    d[i + 1] = g + p;
                    g = c * r - b;

                    // ── accumulate Givens rotation into Q ──
                    for (let k = 0; k < n; k++) {
                        const t      = Q[k][i + 1];
                        Q[k][i + 1]  = s * Q[k][i] + c * t;
                        Q[k][i]      = c * Q[k][i] - s * t;
                    }
                }

                if (!underflow) {
                    d[l] -= p;
                    e[l]  = g;
                    e[m]  = 0.0;
                }
            }
        }

        /* ================================================================
        6.  Sort eigenvalues ascending, pack nEig smallest pairs
        ================================================================ */
        const idx = Array.from({ length: n }, (_, i) => i);
        idx.sort((a, b) => d[a] - d[b]);

        const eigenvalues  = new Array(nEig);
        const eigenvectors = new Array(nEig);

        for (let j = 0; j < nEig; j++) {
            const col = idx[j];
            eigenvalues[j] = d[col];

            const vec = new Array(n);
            for (let k = 0; k < n; k++) vec[k] = Q[k][col];
            eigenvectors[j] = vec;
        }

        return { eigenvalues, eigenvectors };
    }    

    function Eig_Real_General(A, nEig, tol) {
        const n = A.length;
        
        const eps = (tol != null) ? tol : 1e-8;

        // CHANGE 8 (v4): Deflation threshold at machine-epsilon scale, separate
        // from user-facing `eps`. Used for STRUCTURAL zero tests only.
        const epsDefl = Math.max(n * Number.EPSILON, 1e-14);

        // CHANGE 5 (v3): Removed the per-block iteration cap ...
        const maxIterGlobal = 100 * n;


        // LAPACK dlahqr constants used by the exceptional-shift formulas below.
        const dat1 =  0.75;
        const dat2 = -0.4375;

        // 1. Hessenberg reduction (unchanged)
        let H = A.map(row => row.slice());
        let Q = Array.from({length: n}, () => Array(n).fill(0));
        for (let i = 0; i < n; i++) Q[i][i] = 1.0;
        for (let k = 0; k < n - 2; k++) {
            const m = n - k - 1;
            const x = H.slice(k+1, k+1+m).map(r => r[k]);
            let xNorm = Math.sqrt(x.reduce((s,v)=>s+v*v,0));
            if (xNorm < eps) continue;
            const alpha = x[0] >= 0 ? -xNorm : xNorm;
            let v = x.slice();
            v[0] -= alpha;
            let vNorm = Math.sqrt(v.reduce((s,vi)=>s+vi*vi,0));
            if (vNorm < eps) continue;
            v = v.map(vi => vi / vNorm);
            for (let j = 0; j < n; j++) {
                let dot = v.reduce((s,vi,i) => s + vi * H[k+1+i][j], 0) * 2;
                for (let i = 0; i < m; i++) H[k+1+i][j] -= dot * v[i];
            }
            for (let i = 0; i < n; i++) {
                let dot = v.reduce((s,vj,j) => s + H[i][k+1+j] * vj, 0) * 2;
                for (let j = 0; j < m; j++) H[i][k+1+j] -= dot * v[j];
            }
            for (let i = 0; i < n; i++) {
                let dot = v.reduce((s,vj,j) => s + Q[i][k+1+j] * vj, 0) * 2;
                for (let j = 0; j < m; j++) Q[i][k+1+j] -= dot * v[j];
            }
        }

        // 2. Francis double-shift QR iterations with LAPACK-style exceptional shifts
        let hi = n - 1;
        let totalIter = 0;

        // CHANGE 4 (v3): Track iterations spent at the current hi. Reset every
        // time hi advances (i.e., a deflation succeeds). This counter drives
        // exceptional shift injection below.
        let iterAtHi = 0;
        let prevHi   = -1;

        while (hi > 0 && totalIter < maxIterGlobal) {
            if (hi !== prevHi) {
                iterAtHi = 0;
                prevHi   = hi;
            }

            
            let lo = hi;
            // CHANGE 9 (v4): Use epsDefl (structural) not eps (user tol)
            while (lo > 0 && Math.abs(H[lo][lo-1]) > epsDefl * (Math.abs(H[lo-1][lo-1]) + Math.abs(H[lo][lo]))) {
                lo--;
            }

            if (lo === hi) {
                hi--;
                continue;
            }
            if (lo === hi - 1) {
                // Only treat this 2x2 block as converged/terminal if its
                // eigenvalues are genuinely complex (irreducible over the
                // reals). If they are real, deflate immediately with a
                // closed-form rotation below -- do NOT fall through to the
                // shift step, since when lo === hi-1 the shift s,t is built
                // from this same 2x2 block's own trace/determinant, which
                // makes x and y come out to be exactly 0 algebraically. That
                // leaves norm < eps, which skips the update entirely, so H
                // never changes and the loop stalls on this block forever --
                // burning the whole iteration budget and leaving everything
                // above it un-reduced, which corrupts eigenvalues/eigenvectors
                // extracted from those rows afterward.
                const trB   = H[lo][lo] + H[hi][hi];
                const dtB   = H[lo][lo] * H[hi][hi] - H[lo][hi] * H[hi][lo];
                const discB = trB * trB - 4 * dtB;
                if (discB < 0) {
                    hi -= 2;
                    continue;
                }
                // Two distinct REAL eigenvalues: deflate directly with a
                // closed-form similarity rotation R = [v, v_perp], where v is
                // the (real) eigenvector of lambda1 within the block -- the
                // same approach used later for eigenvector back-substitution.
                {
                    const s2   = Math.sqrt(discB) / 2;
                    const lam1 = trB / 2 + s2;
                    const a = H[lo][lo], b = H[lo][hi], c = H[hi][lo], d = H[hi][hi];
                    let vx, vy;
                    if (Math.abs(b) >= Math.abs(c)) { vx = b;           vy = lam1 - a; }
                    else                             { vx = lam1 - d;   vy = c;        }
                    const vnorm = Math.sqrt(vx * vx + vy * vy);
                    if (vnorm > eps) {
                        vx /= vnorm; vy /= vnorm;
                        for (let col = 0; col < n; col++) {
                            const h0 = H[lo][col], h1 = H[hi][col];
                            H[lo][col] =  vx * h0 + vy * h1;
                            H[hi][col] = -vy * h0 + vx * h1;
                        }
                        for (let row = 0; row < n; row++) {
                            const h0 = H[row][lo], h1 = H[row][hi];
                            H[row][lo] = h0 * vx + h1 * vy;
                            H[row][hi] = -h0 * vy + h1 * vx;
                        }
                        for (let row = 0; row < n; row++) {
                            const q0 = Q[row][lo], q1 = Q[row][hi];
                            Q[row][lo] = q0 * vx + q1 * vy;
                            Q[row][hi] = -q0 * vy + q1 * vx;
                        }
                    }
                    // CHANGE 3 (v1): Hard-zero the subdiagonal entry to eliminate
                    // residual floating-point noise from the similarity rotation.
                    H[hi][lo] = 0;
                    hi -= 2;
                    continue;
                }
            }

            // CHANGE 6 (v3): LAPACK dlahqr-style exceptional shifts at iterations
            // 10 and 20 within the same block. The two shifts are anchored at
            // opposite ends of the active window (hi vs. lo), which breaks
            // different stall patterns. A third random-jitter shift at
            // iterations 30, 40, 50, ... handles the very rare cases where
            // both LAPACK shifts also stall.
            let s, t;
            if (iterAtHi === 10) {
                // Exceptional shift 1: anchored at the bottom of the active window.
                const s_ex   = Math.abs(H[hi][hi-1]) + (hi > 1 ? Math.abs(H[hi-1][hi-2]) : 0);
                const hShift = dat1 * s_ex + H[hi][hi];
                s = 2 * hShift;
                t = hShift * hShift - dat2 * s_ex * s_ex;
            } else if (iterAtHi === 20) {
                // Exceptional shift 2: anchored at the top of the active window.
                const s_ex   = Math.abs(H[lo+1][lo]) + (lo+2 <= hi ? Math.abs(H[lo+2][lo+1]) : 0);
                const hShift = dat1 * s_ex + H[lo][lo];
                s = 2 * hShift;
                t = hShift * hShift - dat2 * s_ex * s_ex;
            } else if (iterAtHi > 30 && iterAtHi % 10 === 0) {
                // Exceptional shift 3: random-jitter fallback. Very rare, but
                // guarantees eventual progress on truly pathological matrices
                // where both LAPACK exceptional shifts also fail to reduce.
                const scale  = Math.abs(H[hi][hi]) + Math.abs(H[lo][lo]) + 1;
                const jitter = scale * (Math.random() - 0.5);
                const hShift = H[hi][hi] + jitter;
                s = 2 * hShift;
                t = hShift * hShift + 0.1 * scale * scale;
            } else {
                // Standard Wilkinson double shift from bottom 2x2
                s = H[hi-1][hi-1] + H[hi][hi];
                t = H[hi-1][hi-1] * H[hi][hi] - H[hi-1][hi] * H[hi][hi-1];
            }

            let x = H[lo][lo]*H[lo][lo] + H[lo][lo+1]*H[lo+1][lo] - s*H[lo][lo] + t;
            let y = H[lo+1][lo] * (H[lo][lo] + H[lo+1][lo+1] - s);
            let z = (lo+2 <= hi) ? H[lo+2][lo+1]*H[lo+1][lo] : 0;

            for (let k = lo; k < hi; k++) {
                const vLen = (k < hi - 1) ? 3 : 2;
                let v = [x, y, z].slice(0, vLen);
                let norm = Math.sqrt(v.reduce((a,b)=>a+b*b,0));
                if (norm < eps) {
                    if (k < hi-1) {
                        x = H[k+1][k];
                        y = H[k+2][k];
                        z = (k+3<=hi) ? H[k+3][k] : 0;
                    }
                    continue;
                }
                const alpha = v[0] >= 0 ? -norm : norm;
                v[0] -= alpha;
                norm = Math.sqrt(v.reduce((a,b)=>a+b*b,0));
                v = v.map(vi => vi / norm);
                const r0 = Math.max(k, lo);
                // Left multiply
                for (let j = 0; j < n; j++) {
                    let dot = v.reduce((s,vi,i) => s + vi * H[r0+i][j], 0) * 2;
                    for (let i = 0; i < vLen; i++) H[r0 + i][j] -= dot * v[i];
                }
                // Right multiply
                for (let i = 0; i < n; i++) {
                    let dot = v.reduce((s,vj,j) => s + H[i][k+j] * vj, 0) * 2;
                    for (let j = 0; j < vLen; j++) H[i][k + j] -= dot * v[j];
                }
                // Accumulate Q
                for (let i = 0; i < n; i++) {
                    let dot = v.reduce((s,vj,j) => s + Q[i][k+j] * vj, 0) * 2;
                    for (let j = 0; j < vLen; j++) Q[i][k + j] -= dot * v[j];
                }
                if (k < hi - 1) {
                    x = H[k + 1][k];
                    y = H[k + 2][k];
                    z = (k + 3 <= hi) ? H[k + 3][k] : 0;
                }
            }
            totalIter++;
            iterAtHi++;
        }

        // Optional diagnostic: warn if the global budget was exhausted. On any
        // reasonable random matrix this should not fire. If it does, save the
        // matrix off -- it is genuinely pathological.
        if (totalIter >= maxIterGlobal && hi > 0) {
            console.warn(`Eig_Real_General: iteration budget exhausted with hi=${hi}. ` +
                        `Eigenvalues near this position may be inaccurate.`);
        }

        // 3. Extract eigenvalues
        // CHANGE 7 (v2): Use a RELATIVE subdiagonal threshold that matches the
        // criterion used by the QR iteration itself when searching for `lo`.
        // The old check `< eps` was absolute, which meant a subdiagonal that had
        // been deflated during iteration by the relative criterion (e.g. 1e-5
        // when diagonals are ~1000) could still look "not negligible" to the
        // extraction phase, spuriously triggering the 2x2 block path on what is
        // really an already-deflated 1x1 eigenvalue.
        const eigenvalues = [];
        const schurIdx = [];
        let i = 0;
        
        while (i < n) {
            // CHANGE 10 (v4): Use epsDefl (structural) not eps (user tol) -- must
            // match the iteration's criterion in section 2.
            const subdiagNegligible = (i === n - 1) ||
                (Math.abs(H[i + 1][i]) < epsDefl * (Math.abs(H[i][i]) + Math.abs(H[i + 1][i + 1])));

            if (subdiagNegligible) {
                eigenvalues.push(H[i][i]);
                schurIdx.push({type: 'r', idx: i});
                i++;
            } else {
                const tr = H[i][i] + H[i + 1][i + 1];
                const dt = H[i][i] * H[i + 1][i + 1] - H[i][i + 1] * H[i + 1][i];
                const re = tr / 2;
                const disc = tr * tr - 4 * dt;
                if (disc >= 0) {
                    // Undeflated 2x2 block with two DISTINCT REAL eigenvalues.
                    // Should be rare now that the main loop above deflates real
                    // 2x2 blocks directly, but kept as a safety net -- deflate
                    // this block directly with a single closed-form similarity
                    // rotation R = [v, v_perp], where v is the (real)
                    // eigenvector of lambda1 within the block. This exactly
                    // zeros H[i+1][i] and also cleans up the coupling that later
                    // blocks' eigenvector solves read through.
                    const s = Math.sqrt(disc) / 2;
                    const lambda1 = re + s, lambda2 = re - s;
                    const a = H[i][i], b = H[i][i + 1], c = H[i + 1][i], d = H[i + 1][i + 1];
                    let vx, vy;
                    if (Math.abs(b) >= Math.abs(c)) { vx = b;            vy = lambda1 - a; }
                    else                             { vx = lambda1 - d; vy = c;            }
                    const vnorm = Math.sqrt(vx * vx + vy * vy);
                    if (vnorm > eps) {
                        vx /= vnorm; vy /= vnorm;
                        for (let col = 0; col < n; col++) {
                            const h0 = H[i][col], h1 = H[i + 1][col];
                            H[i][col]     =  vx * h0 + vy * h1;
                            H[i + 1][col] = -vy * h0 + vx * h1;
                        }
                        for (let row = 0; row < n; row++) {
                            const h0 = H[row][i], h1 = H[row][i + 1];
                            H[row][i]     = h0 * vx + h1 * vy;
                            H[row][i + 1] = -h0 * vy + h1 * vx;
                        }
                        for (let row = 0; row < n; row++) {
                            const q0 = Q[row][i], q1 = Q[row][i + 1];
                            Q[row][i]     = q0 * vx + q1 * vy;
                            Q[row][i + 1] = -q0 * vy + q1 * vx;
                        }
                    }
                    // CHANGE 2 (v1): Hard-zero the subdiagonal entry to eliminate
                    // residual floating-point noise from the similarity rotation.
                    H[i + 1][i] = 0;
                    eigenvalues.push(lambda1);
                    eigenvalues.push(lambda2);
                    // CHANGE 1 (v1): Route both real eigenvalues through the REAL
                    // back-substitution branch, and give each its own seed
                    // position within the (now triangular) block. Previously
                    // both were pushed as {type: 'c', idx: i}, which sent them
                    // through the complex branch and caused lambda2 to receive
                    // a numerically unstable seed for a degenerate local 2x2
                    // solve -- producing catastrophic residual failure on pairs
                    // of close real eigenvalues.
                    schurIdx.push({type: 'r', idx: i});      // lambda1 seeds at row i
                    schurIdx.push({type: 'r', idx: i + 1});  // lambda2 seeds at row i+1
                } else {
                    // Genuine complex-conjugate pair
                    const im = Math.sqrt(-disc) / 2;
                    eigenvalues.push(new ComplexNum(re, im));
                    eigenvalues.push(new ComplexNum(re, -im));
                    schurIdx.push({type: 'c', idx: i});
                    schurIdx.push({type: 'c', idx: i});
                }
                i += 2;
            }
        }

        // 4. Back-substitution for eigenvectors
        const eigenvectors = [];
        for (let j = 0; j < eigenvalues.length; j++) {
            const si = schurIdx[j];
            if (si.type === 'r') {
                const lambda = eigenvalues[j];
                let y = new Array(n).fill(0); y[si.idx] = 1;
                let r = si.idx - 1;
                while (r >= 0) {
                    // If rows (r-1, r) form a genuinely undeflated 2x2 block
                    // (a complex-conjugate pair that can't be triangularized
                    // over the reals), H[r][r-1] is NOT negligible. Solving row
                    // r alone as a simple scalar divide silently drops the
                    // H[r][r-1]*y[r-1] coupling term, giving a wrong y[r] --
                    // and since y[r-1] is computed next FROM that wrong y[r],
                    // the error propagates through the rest of the vector.
                    // Solve both rows of the block together instead.
                    if (r > 0 && Math.abs(H[r][r-1]) > eps * (Math.abs(H[r-1][r-1]) + Math.abs(H[r][r]))) {
                        let sum0 = 0, sum1 = 0;
                        for (let c = r + 1; c <= si.idx; c++) {
                            sum0 += H[r-1][c] * y[c];
                            sum1 += H[r][c] * y[c];
                        }
                        const a = H[r-1][r-1] - lambda, b = H[r-1][r];
                        const cc = H[r][r-1],            d = H[r][r] - lambda;
                        const det = a * d - b * cc;
                        if (Math.abs(det) > eps) {
                            y[r-1] = (-sum0 * d + b * sum1) / det;
                            y[r]   = (-a * sum1 + cc * sum0) / det;
                        } else {
                            y[r-1] = 0; y[r] = 0;
                        }
                        r -= 2;
                    } else {
                        let sum = 0;
                        for (let c = r + 1; c <= si.idx; c++) sum += H[r][c] * y[c];
                        const d = H[r][r] - lambda;
                        y[r] = Math.abs(d) > eps ? -sum / d : 0;
                        r -= 1;
                    }
                }
                let v = new Array(n).fill(0);
                for (let r = 0; r < n; r++) for (let c = 0; c <= si.idx; c++) v[r] += Q[r][c] * y[c];
                const norm = Math.sqrt(v.reduce((a,b)=>a+b*b,0) || 1);
                v = v.map(x => x / norm);
                eigenvectors.push(v);
            } else {
                // 2x2 undeflated block: solve independently for THIS eigenvalue
                // (works whether lambda is real or a genuine complex-conjugate
                // partner -- do not assume the second eigenvector is just the
                // conjugate of the first, since for a real-eigenvalue block the
                // two eigenvalues are NOT conjugates of one another).
                const blk = si.idx;
                const lambda = eigenvalues[j];
                let y = new Array(n).fill(null).map(() => new ComplexNum(0, 0));
                // Solve the (singular) 2x2 block (H_block - lambda*I) y = 0.
                // Row blk:   c00*y[blk] + c01*y[blk+1] = 0
                // Row blk+1: c10*y[blk] + c11*y[blk+1] = 0
                // The two rows are proportional (lambda is an eigenvalue of the
                // block); pick whichever row/pivot has the larger magnitude to
                // divide by for numerical stability. This matters in particular
                // once a block has been fully deflated to triangular: for its
                // TOP eigenvalue, c00 is exactly 0, so fixing y[blk+1]=1 and
                // solving for y[blk] via c00 would divide by zero -- here we
                // instead recognize that case and set y[blk]=1, y[blk+1]=0.
                let c00 = new ComplexNum(H[blk][blk], 0).Subtract(lambda);
                let c01 = new ComplexNum(H[blk][blk + 1], 0);
                let c10 = new ComplexNum(H[blk + 1][blk], 0);
                let c11 = new ComplexNum(H[blk + 1][blk + 1], 0).Subtract(lambda);
                if (c01.Abs() >= c00.Abs() && c01.Abs() > eps) {
                    y[blk] = new ComplexNum(1, 0);
                    y[blk + 1] = c00.Multiply(-1).Divide(c01);
                } else if (c00.Abs() > eps) {
                    y[blk + 1] = new ComplexNum(1, 0);
                    y[blk] = c01.Multiply(-1).Divide(c00);
                } else if (c11.Abs() >= c10.Abs() && c11.Abs() > eps) {
                    y[blk] = new ComplexNum(1, 0);
                    y[blk + 1] = c10.Multiply(-1).Divide(c11);
                } else if (c10.Abs() > eps) {
                    y[blk + 1] = new ComplexNum(1, 0);
                    y[blk] = c11.Multiply(-1).Divide(c10);
                } else {
                    // Fully diagonal block: this eigenvalue occupies row blk
                    y[blk] = new ComplexNum(1, 0);
                }
                let r = blk - 1;
                while (r >= 0) {
                    // If rows (r-1, r) form a genuinely undeflated 2x2 block
                    // (another complex-conjugate pair that couldn't be
                    // triangularized over the reals -- this happens whenever a
                    // matrix has more than one non-real eigenvalue pair, since
                    // back-substitution for the topmost pair's eigenvector must
                    // pass through any lower pairs' still-undeflated blocks),
                    // H[r][r-1] is NOT negligible. Solving row r alone as a
                    // simple scalar divide silently drops the H[r][r-1]*y[r-1]
                    // coupling term, giving a wrong y[r] -- and since y[r-1] is
                    // computed next FROM that wrong y[r], the error propagates
                    // through the rest of the vector. Solve both rows of the
                    // block together instead (same approach as the real-
                    // eigenvalue branch above, done here in complex arithmetic
                    // since lambda is complex).
                    if (r > 0 && Math.abs(H[r][r - 1]) > eps * (Math.abs(H[r - 1][r - 1]) + Math.abs(H[r][r]))) {
                        let sum0 = new ComplexNum(0, 0), sum1 = new ComplexNum(0, 0);
                        for (let c = r + 1; c <= blk + 1; c++) {
                            sum0 = sum0.Add(new ComplexNum(H[r - 1][c], 0).Multiply(y[c]));
                            sum1 = sum1.Add(new ComplexNum(H[r][c], 0).Multiply(y[c]));
                        }
                        const a  = new ComplexNum(H[r - 1][r - 1], 0).Subtract(lambda);
                        const b  = new ComplexNum(H[r - 1][r], 0);
                        const cc = new ComplexNum(H[r][r - 1], 0);
                        const d  = new ComplexNum(H[r][r], 0).Subtract(lambda);
                        const det = a.Multiply(d).Subtract(b.Multiply(cc));
                        if (det.Abs() > eps) {
                            y[r - 1] = sum0.Multiply(-1).Multiply(d).Add(b.Multiply(sum1)).Divide(det);
                            y[r]     = a.Multiply(sum1.Multiply(-1)).Add(cc.Multiply(sum0)).Divide(det);
                        } else {
                            y[r - 1] = new ComplexNum(0, 0);
                            y[r]     = new ComplexNum(0, 0);
                        }
                        r -= 2;
                    } else {
                        let sum = new ComplexNum(0, 0);
                        for (let c = r + 1; c <= blk + 1; c++) sum = sum.Add(new ComplexNum(H[r][c], 0).Multiply(y[c]));
                        const d = new ComplexNum(H[r][r], 0).Subtract(lambda);
                        y[r] = d.Abs() > eps ? sum.Multiply(-1).Divide(d) : new ComplexNum(0, 0);
                        r -= 1;
                    }
                }
                let v = new Array(n).fill(null).map(() => new ComplexNum(0, 0));
                for (let r = 0; r < n; r++) {
                    for (let c = 0; c <= blk + 1; c++) v[r] = v[r].Add(y[c].Multiply(Q[r][c]));
                }
                let nrm = Math.sqrt(v.reduce((s, z) => s + z.Re*z.Re + z.Im*z.Im, 0) || 1);
                v = v.map(z => z.Multiply(1 / nrm));
                // Fix the arbitrary complex phase: rotate so the largest-magnitude
                // component is real and positive (matches MATLAB/LAPACK convention).
                // Complex eigenvectors are only unique up to a unit-magnitude phase
                // factor e^{i*theta}; without this, the vector returned here is a
                // valid eigenvector but its phase depends on arbitrary numerical
                // choices made during the back-substitution above.
                {
                    let maxAbs = -1, maxIdx = 0;
                    for (let r = 0; r < n; r++) {
                        const m = v[r].Abs();
                        if (m > maxAbs) { maxAbs = m; maxIdx = r; }
                    }
                    if (maxAbs > eps) {
                        const theta = v[maxIdx].Angle();
                        const rot = new ComplexNum(Math.cos(-theta), Math.sin(-theta));
                        v = v.map(z => z.Multiply(rot));
                    }
                }
                eigenvectors.push(v);
            }
        }

        // 5. Return in natural Schur order (matches MATLAB/LAPACK dgeev ordering).
        // MATLAB/LAPACK do NOT sort eigenvalues by magnitude -- they return them
        // in the order they fall out of the real Schur form produced by the
        // Francis double-shift QR algorithm. Since this function uses the same
        // algorithm, the order eigenvalues/eigenvectors were pushed into the
        // arrays above already matches that convention. Sorting by |lambda|
        // here was reordering away from that convention -- removed.
        const nRet = Math.min(nEig, eigenvalues.length);
        const retVals = eigenvalues.slice(0, nRet);
        const retVecs = eigenvectors.slice(0, nRet);
        return { eigenvalues: retVals, eigenvectors: retVecs };
    }

    function Eig_Complex_Hermitian_Tridiag(A, nEig, tol)  { return null; }
    function Eig_Complex_Hermitian_PosDef(A, nEig, tol)   { return null; }
    function Eig_Complex_Hermitian(A, nEig, tol)          { return null; }
    function Eig_Complex_General(A, nEig, tol)            { return null; }

}
//-------------------------------------------------------------------------------------------------


// let matrix = [
//     [0.0459957881361126, 0.9319231629176088, -0.6731849382095623, -0.06631617676504176, -0.054151321164552346, 0.5381041137226021, 0.14031534836598247],
//     [0.7230806317786072, 0.8389739684959197, -0.5559112709216263, 0.2907508236706531, 0.5389777352804297, -0.30787960315479346, -0.024921545828543223],
//     [-0.4150705692917882, 0.2579891861170427, 0.8491393731905861, 0.564270850325866, 0.18740694508631117, 0.6437072729531021, -0.500340024253247],
//     [0.11379715707422733, -0.12483176462884349, -0.935497349888254, 0.2034144422117179, 0.7888539910253893, 0.8676234091425021, 0.9789812981861912],
//     [0.09087789721843542, -0.0008440116281205334, -0.5590170062503372, -0.12500838148917737, -0.7610945170980536, 0.5091476901691141, 0.2735707764859763],
//     [0.7563045015787253, -0.5524659086965296, -0.9691530721863604, 0.019788332404966802, 0.31398173551410835, -0.23562711133527547, 0.16894912344341084],
//     [0.34054481209664433, -0.10808983158788732, -0.6975933610157565, -0.05603176974972013, 0.25601665571613275, 0.21831334753936216, 0.3895719735263645],
// ];


// let res = Eig(matrix, { auto: false, isReal: true, isSymmHermi: false, isPosDef: false, isTridiagonal: false, tol: 1e-14 });
// console.log(res.eigenvalues)
// console.log(res.eigenvectors)

// let ver = Verify_Eig(matrix, res.eigenvalues, res.eigenvectors, 1e-14)
// console.log(ver)

// console.log(ver.residual);
// console.log(ver.perEigenPass);

//-------------------------------------------------------------------------------------------------




