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
    let tol = 1e-6;
    let NumSim = 500;
    let mMax = 5, mMin=100;

    // Counters
    let nPass = 0, nFail = 0;
    const failures = [];

    for (let rr=0; rr<NumSim; rr++) {

        // Randomly select the size of matrix
        let mSize = Math.floor(Math.random() * (mMax - mMin + 1)) + mMin;

        // Randomly choose real or complex and Build a random base matrix-B
        isReal         = randomBoolean();
        B              = isReal ? Create_Real_Matrix(mSize) : Create_Complex_Matrix(mSize);

        // Randomly choose matrix class; construction follows the flag
        isSymmHermi    = randomBoolean();
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
        if (res == null || res.eigenvalues == null || res.eigenvectors == null) {
            continue;
        }
        
        // ---- Verification ----------------------------------------------
        //{ pass, residual, ortho, details }
        let Ver_result = Verify_Eig(A, res.eigenvalues, res.eigenvectors, tol)

        console.log(rr)
        console.log('IsReal               : ' + isReal)
        console.log('isSymmHermi          : ' + isSymmHermi)
        console.log('isPosDef             : ' + isPosDef)
        console.log('isTridiagonal        : ' + isTridiagonal)
        console.log('Ver_result.pass      : ' + Ver_result.pass)
        console.log('Ver_result.residual  : ' + Ver_result.pass)
        console.log('Ver_result.ortho     : ' + Ver_result.ortho)
        console.log('Ver_result.details   : ' + Ver_result.details)

        if (Ver_result.pass == false) { return 'Error'; }

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


    function Verify_Eig(A, eigenvalues, eigenvectors, tol) {

        if (tol == null) tol = 1e-6;

        const nEig    = eigenvalues.length;
        const n       = A.length;
        const details = [];
        let   pass    = true;

        // ── 1. Build V (n × nEig) from the eigenvectors array ──
        //    eigenvectors[j] is the j-th eigenvector of length n
        //    V[i][j] = eigenvectors[j][i]   (columns = eigenvectors)
        const V = [];
        for (let i = 0; i < n; i++) {
            V[i] = [];
            for (let j = 0; j < nEig; j++) {
                V[i][j] = eigenvectors[j][i];
            }
        }

        // ── 2. Residual:  ‖A·V − V·diag(λ)‖_F  /  ‖A‖_F ──
        //    A·V  is n × nEig
        //    V·diag(λ) means: column j of V scaled by λ[j]
        const AV = Multiply(A, V);

        const VL = [];
        for (let i = 0; i < n; i++) {
            VL[i] = [];
            for (let j = 0; j < nEig; j++) {
                VL[i][j] = Multiply(V[i][j], eigenvalues[j]);
            }
        }


        const Res_mat  = Subtract(AV, VL);
        const normA    = Norm(A, 'fro');
        const normRes  = Norm(Res_mat, 'fro');
        const residual = normRes / (normA > 0 ? normA : 1);

        if (residual >= tol) {
            pass = false;
            details.push('FAIL residual = ' + residual.toExponential(3) + ' >= ' + tol);
        } else {
            details.push('OK   residual = ' + residual.toExponential(3));
        }

        // ── 3. Orthonormality:  ‖Vᵀ·V − I‖_F ──
        const VtV     = Multiply(Transpose(V), V);           // nEig × nEig
        const I_nEig  = Eye(nEig);
        const Ort_mat = Subtract(VtV, I_nEig);
        const ortho   = Norm(Ort_mat, 'fro');

        if (ortho >= tol) {
            pass = false;
            details.push('FAIL ortho    = ' + ortho.toExponential(3) + ' >= ' + tol);
        } else {
            details.push('OK   ortho    = ' + ortho.toExponential(3));
        }

        // ── 4. Individual eigenvector unit-norm check ──
        for (let j = 0; j < nEig; j++) {
            const vj_norm = Norm(eigenvectors[j], 2);
            if (Math.abs(vj_norm - 1.0) >= tol) {
                pass = false;
                details.push('FAIL eigenvector[' + j + '] norm = ' + vj_norm.toExponential(3));
            }
        }

        // ── 5. Eigenvalue ordering (ascending) ──
        for (let j = 1; j < nEig; j++) {
            if (eigenvalues[j] < eigenvalues[j - 1] - tol) {
                pass = false;
                details.push('FAIL eigenvalues not sorted: λ[' + (j-1) + ']='
                    + eigenvalues[j-1].toExponential(3) + ' > λ[' + j + ']='
                    + eigenvalues[j].toExponential(3));
            }
        }

        // ── 6. Trace check (only when all eigenvalues extracted) ──
        if (nEig === n) {
            let traceA  = 0;
            let sumEig  = 0;
            for (let i = 0; i < n; i++) { traceA += A[i][i]; }
            for (let j = 0; j < nEig; j++) { sumEig += eigenvalues[j]; }

            const traceErr = Math.abs(traceA - sumEig) / (Math.abs(traceA) > 0 ? Math.abs(traceA) : 1);
            if (traceErr >= tol) {
                pass = false;
                details.push('FAIL trace mismatch: tr(A)=' + traceA.toExponential(3)
                    + '  Σλ=' + sumEig.toExponential(3));
            } else {
                details.push('OK   trace match   = ' + traceErr.toExponential(3));
            }
        }

        return { pass, residual, ortho, details };
    }
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
    //   Eig([[2, 1], [1, 2]])
    //       → { Val: [1, 3], Vec: [[-0.707, 0.707], [0.707, 0.707]] }
    //   Eig([[4, 1, 0], [1, 3, 1], [0, 1, 2]], { nEig: 2 })
    //       → { Val: [λ₁, λ₂], Vec: [[...], [...]] }  — 2 smallest-magnitude eigenpairs
    //   Eig(A, { auto: false, isReal: true, isSymmHermi: true, isPosDef: false,
    //            isTridiagonal: false, nEig: 3 })
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
        else                                   { return Eig_Real_General(A, nEigResolved);       }  // Real general

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





function Eig_Real_General(A, nEig) {
    const n = A.length;
    const eps = 1e-8;
    const maxIter = 300 * n;

    // 1. Hessenberg reduction
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

    // 2. Repeated QR iterations with shifts for convergence to Schur form
    let hi = n - 1;
    let totalIter = 0;

    while (hi > 0 && totalIter < maxIter) {
        let lo = hi;
        while (lo > 0 && Math.abs(H[lo][lo-1]) > eps * (Math.abs(H[lo-1][lo-1]) + Math.abs(H[lo][lo]))) {
            lo--;
        }

        if (lo === hi) {
            hi--;
            continue;
        }
        if (lo === hi - 1) {
            hi -= 2;
            continue;
        }

        // Shift (use bottom 2x2)
        const s = H[hi-1][hi-1] + H[hi][hi];
        const t = H[hi-1][hi-1] * H[hi][hi] - H[hi-1][hi] * H[hi][hi-1];

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
    }

    // 3. Extract eigenvalues
    const eigenvalues = [];
    const schurIdx = [];
    let i = 0;
    while (i < n) {
        if (i === n - 1 || Math.abs(H[i + 1][i]) < eps) {
            eigenvalues.push(H[i][i]);
            schurIdx.push({type: 'r', idx: i});
            i++;
        } else {
            const tr = H[i][i] + H[i + 1][i + 1];
            const dt = H[i][i] * H[i + 1][i + 1] - H[i][i + 1] * H[i + 1][i];
            const re = tr / 2;
            const im = 0.5 * Math.sqrt(Math.max(0, tr*tr - 4*dt));
            eigenvalues.push(new ComplexNum(re, im));
            eigenvalues.push(new ComplexNum(re, -im));
            schurIdx.push({type: 'c', idx: i});
            schurIdx.push({type: 'c', idx: i});
            i += 2;
        }
    }

    // 4. Back-substitution for eigenvectors
    const eigenvectors = [];
    const done = new Set();

    for (let j = 0; j < eigenvalues.length; j++) {
        const si = schurIdx[j];
        if (si.type === 'r') {
            const lambda = eigenvalues[j];
            let y = new Array(n).fill(0); y[si.idx] = 1;
            for (let r = si.idx - 1; r >= 0; r--) {
                let sum = 0;
                for (let c = r + 1; c <= si.idx; c++) sum += H[r][c] * y[c];
                const d = H[r][r] - lambda;
                y[r] = Math.abs(d) > eps ? -sum / d : 0;
            }
            let v = new Array(n).fill(0);
            for (let r = 0; r < n; r++) for (let c = 0; c <= si.idx; c++) v[r] += Q[r][c] * y[c];
            const norm = Math.sqrt(v.reduce((a,b)=>a+b*b,0) || 1);
            v = v.map(x => x / norm);
            eigenvectors.push(v);
        } else if (!done.has(si.idx)) {
            done.add(si.idx);
            const blk = si.idx;
            const lambda = eigenvalues[j];
            let y = new Array(n).fill(null).map(() => new ComplexNum(0, 0));
            y[blk + 1] = new ComplexNum(1, 0);

            let c00 = new ComplexNum(H[blk][blk], 0).Subtract(lambda);
            let c01 = new ComplexNum(H[blk][blk + 1], 0);
            if (c00.Abs() > eps) y[blk] = c01.Multiply(-1).Divide(c00);

            for (let r = blk - 1; r >= 0; r--) {
                let sum = new ComplexNum(0, 0);
                for (let c = r + 1; c <= blk + 1; c++) sum = sum.Add(new ComplexNum(H[r][c], 0).Multiply(y[c]));
                let d = new ComplexNum(H[r][r], 0).Subtract(lambda);
                y[r] = d.Abs() > eps ? sum.Multiply(-1).Divide(d) : new ComplexNum(0, 0);
            }

            let v = new Array(n).fill(null).map(() => new ComplexNum(0, 0));
            for (let r = 0; r < n; r++) {
                for (let c = 0; c <= blk + 1; c++) v[r] = v[r].Add(y[c].Multiply(Q[r][c]));
            }
            let nrm = Math.sqrt(v.reduce((s, z) => s + z.Re*z.Re + z.Im*z.Im, 0) || 1);
            v = v.map(z => z.Multiply(1 / nrm));

            eigenvectors.push(v);
            eigenvectors.push(v.map(z => z.Conj()));
        }
    }

    // 5. Sort by |λ| and take smallest nEig
    const idx = Array.from({length: eigenvalues.length}, (_,k)=>k)
        .sort((a,b) => {
            const ma = eigenvalues[a] instanceof ComplexNum ? eigenvalues[a].Abs() : Math.abs(eigenvalues[a]);
            const mb = eigenvalues[b] instanceof ComplexNum ? eigenvalues[b].Abs() : Math.abs(eigenvalues[b]);
            return ma - mb;
        });

    const retVals = [], retVecs = [];
    for (let j = 0; j < Math.min(nEig, eigenvalues.length); j++) {
        retVals.push(eigenvalues[idx[j]]);
        retVecs.push(eigenvectors[idx[j]]);
    }

    return { eigenvalues: retVals, eigenvectors: retVecs };
}




    function Eig_Complex_Hermitian_Tridiag(A, nEig, tol)  { return null; }
    function Eig_Complex_Hermitian_PosDef(A, nEig, tol)   { return null; }
    function Eig_Complex_Hermitian(A, nEig, tol)          { return null; }
    function Eig_Complex_General(A, nEig, tol)            { return null; }
}






// Eig(A)
// │
// ├── Scalar / ComplexNum / 1×1        → closed form (trivial)
// ├── 2×2                              → closed form (quadratic formula)
// │
// ├── isDiagonal                       → O(n) read diagonal
// ├── isTriangular                     → O(n) diagonal + back-substitution
// │
// ├── isReal
// │   ├── isSymmHermi
// │   │   ├── isTridiagonal            → [Branch 1a] skip reduction, QR iteration only
// │   │   └── (general symmetric)      → [Branch 1b] Householder → tridiagonal → same QR kernel as 1a
// │   └── (general)                    → [Branch 2]  Householder → Hessenberg → Francis double-shift
// │
// └── isComplex
//     ├── isSymmHermi
//     │   ├── isTridiagonal            → [Branch 3a] phase-absorb → same real QR kernel as 1a
//     │   └── (general Hermitian)      → [Branch 3b] Householder + phase → same real QR kernel as 1a
//     └── (general)                    → [Branch 4]  Householder → Hessenberg → complex single-shift QR