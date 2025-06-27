// stricter parsing and error handling
"use strict";

//-----------------------------------------------------------------------------------------------------------------------------------------------------------
// Math Library
// The dimensions of the function arguments are not checked and must be consistent.
function infoArray(a) {
    if (Array.isArray(a)) {
       if      (typeof(a[0].length) === 'undefined') { return { 'is1D': true,   'is2D': false } }
       else if (typeof(a[0].length) === 'number')    { return { 'is1D': false,  'is2D': true  } }
    }
    else { return { 'is1D': false,  'is2D': false  } }
}
function Abs(a) {
   let i, IA = infoArray(a);
   if      (isComplexNum(a)) { return a.Abs(); }
   else if (isNumber(a))     { return Math.abs(a); }
   else if (IA.is1D)         { return a.map((v,i) => Abs(v) ); }
   else if (IA.is2D)         { return a.map((v,i) => Abs(v) ); }
}
function Add(a, b) {
   let i, IA = infoArray(a), IB = infoArray(b);
   if      (isComplexNum(a) && isComplexNum(b)) { return a.Add(b); }
   else if (isComplexNum(a) && isNumber(b)    ) { return a.Add(b); }
   else if (isNumber(a)     && isComplexNum(b)) { return b.Add(a); }
   else if (isNumber(a)     && isNumber(b)    ) { return a+b;      }
   else if (IA.is1D) {
       if      (isComplexNum(b) || isNumber(b)) { return a.map((v, i) => Add(v, b   ) ); }
       else if (IB.is1D)                        { return a.map((v, i) => Add(v, b[i]) ); }
   }
   else if (IA.is2D) {
       if      (isComplexNum(b) || isNumber(b)) { return a.map((v,i) => Add(v, b   ) ); }
       else if (IB.is2D)                        { return a.map((v,i) => Add(v, b[i]) ); }
   }
}
function Angle(a) {
   let i, IA = infoArray(a);
   if      (isComplexNum(a)) { return a.Angle(); }
   else if (isNumber(a))     { if (a<0) { return Math.PI; } else { return 0; } }
   else if (IA.is1D)         { return a.map((v,i) => Angle(v) ); }
   else if (IA.is2D)         { return a.map((v,i) => Angle(v) ); }
}
function Balance(a, RADIX) {
   // Returns a balanced 2D-Array square matrix with the same eigenvalues as 2D-Array-a[].
   // The parameter RADIX should be the machine’s floating-point radix, and It is taken as 2.0 by default.
   // The idea of balancing is to use similarity transformations to make corresponding rows and columns
   // of the A[,] matrix have comparable norms, thus reducing the overall norm of the matrix while leaving the eigenvalues unchanged.
   // A symmetric matrix is already balanced and is unaffected by this procedure.
   // It is therefore recommended to always balance non-symmetric matrices.

   let last, j, i, n, data=[], s, r, g, f, c, sqrdx;
   let IA = infoArray(a);

   if (RADIX == null) { RADIX=2.0; }

   if      (isNumber(a) || isComplexNum(a)) { return a; }
   else if (IA.is1D || !IA.is2D)            { throw new Error("Matrix must be square."); }
   else if (IA.is2D) {
       // Make a copy of the data[,] matrix
       data  = Copy(a);
       n     = data.length;
       sqrdx = RADIX * RADIX;
       last  = 0;

       while (last == 0) {
           last = 1;
           for (i=0; i<n; i++) { //Calculate row and column norms.
               r = 0.0;
               c = 0.0;
               for (j=0; j<n; j++) {
                   if (j != i) {
                       c += Abs(data[j][i]);
                       r += Abs(data[i][j]);
                   }
               }
               if ((c != 0) && (r != 0)) { // If both are nonzero,
                   g = r / RADIX;
                   f = 1.0;
                   s = c + r;
                   while (c < g) {  // find the integer power of the machine radix that comes closest to balancing the matrix
                       f *= RADIX;
                       c *= sqrdx;
                   }
                   g = r * RADIX;
                   while (c > g) {
                       f /= RADIX;
                       c /= sqrdx;
                   }
                   if (((c + r) / f) < (0.95 * s)) {
                       last = 0;
                       g = 1.0 / f;
                       for (j = 0; j < n; j++) {data[i][j] *= g;}  // Apply similarity transformation
                       for (j = 0; j < n; j++) {data[j][i] *= f;}
                   }
               }
           }
       }
       return data;
   }
}
function CheckRelation(a, rel, th) {
    // Returns true if the relation is correct, and returns false otherwise.
    let IA = infoArray(a), res=[];
    if (isNumber(a)) {
        if      ((rel == '>')  && (a >  th)) { return true;  }
        else if ((rel == '>=') && (a >= th)) { return true;  } 
        else if ((rel == '<')  && (a <  th)) { return true;  } 
        else if ((rel == '<=') && (a <= th)) { return true;  }
        else if ((rel == '==') && (a == th)) { return true;  }
        else if ((rel == '!=') && (a != th)) { return true;  }
        else                                 { return false; }
    }
    else if (isComplexNum(a)) { return CheckRelation(Abs(a), rel, th); }
    else if (IA.is1D) { return a.map((v, i) => CheckRelation(v, rel, th)); }
    else if (IA.is2D) { return a.map((v, i) => CheckRelation(v, rel, th)); }
}
function Chol(a) {
   // A symmetric positive definite matrix is a symmetric matrix with all positive eigenvalues
   // Factorizes symmetric positive definite 2D-Array-a[] into an lower triangular L that satisfies L=L*Transpose(L)
   // 2D-Array-a[] must be square, symmetric, and positive definite.

   let m, n, r, c, sum, j, R=[], IA = infoArray(a);

   // Return if not 2D-Array-a[]
   if (!IA.is2D) { return {'R': [], 'isPosDef': false }; }

   // Size of the matrix A[]
   m = a.length;
   n = a[0].length;

   // Return if not square and symmetric
   if ((n!=m) || !IsSymmetric(a)) { return {'R': [], 'isPosDef': false }; }

   // Choloskey Factorization L[] matrix
   L = Zeros(n, n);

   for (r=0; r<n; r++) {
       for (c=0; c<=r; c++) {
            if (c==r) {
               sum = 0; for (j=0; j<c; j++) { sum += L[c][j] * L[c][j]; }
               if (sum >= a[c][c]) {
                   // Inside the Math.sqrt cannot be zero or negative. if so, it is not positive definite
                   return {'R': [], 'isPosDef': false };
               }
               L[c][c] = Math.sqrt(a[c][c] - sum);
            } else {
               sum = 0;
               for (j=0; j<c; j++) {
                   sum += L[r][j] * L[c][j];
               }
               L[r][c] = (1.0 / L[c][c]) * (a[r][c] - sum);
            }
       }
   }
   return {'R': L, 'isPosDef': true };
}
function Copy(a) {
   let i, IA = infoArray(a);
   if      (isComplexNum(a)) { return new ComplexNum(a.Re, a.Im); }
   if      (isNumber(a))     { return a;                          }
   else if (IA.is1D)         { return a.map((v,i) => Copy(v) );   }
   else if (IA.is2D)         { return a.map((v,i) => Copy(v) );   }
}
function Conj(a) {
   let i, IA = infoArray(a);
   if      (isComplexNum(a)) { return a.Conj(); }
   else if (isNumber(a))     { return a; }
   else if (IA.is1D)         { return a.map((v,i) => Conj(v) ); }
   else if (IA.is2D)         { return a.map((v,i) => Conj(v) ); }
}
function Concat(a, b, dim) {
    // Returns the merge of two arrays
    // If 2D, merges using dim-direction (dim = 0 means Merge by adding new rows,   dim = 1 means Merge by add new columns)

    if (isNumber(a)) {
        if      (IsArrayEmpty(b))    { return [a];         }
        else if (isNumber(b))        { return [a, b];      }
        else if (isComplexNum(b))    { return [a, b];      }
        else if (infoArray(b).is1D)  { return [a, ...b];   }
        else if (infoArray(b).is2D)  { return [[a], ...b]; }
    }
    else if (isComplexNum(a)) {
        if      (IsArrayEmpty(b))    { return [a];       }
        else if (isNumber(b))        { return [a, b];    }
        else if (isComplexNum(b))    { return [a, b];    }
        else if (infoArray(b).is1D)  { return [a, ...b]; }
        else if (infoArray(b).is2D)  { return [[a], ...b]; }
    }
    else if (IsArrayEmpty(a)) {
        if      (IsArrayEmpty(b))    { return [];     }
        else if (isNumber(b))        { return [b];    }
        else if (isComplexNum(b))    { return [b];    }
        else if (infoArray(b).is1D)  { return [...b]; }
        else if (infoArray(b).is2D)  { return [...b]; }
    }
    else if (infoArray(a).is1D) { 
        if      (IsArrayEmpty(b))    { return a;             }
        else if (isNumber(b))        { return [...a, b];     }
        else if (isComplexNum(b))    { return [...a, b];     }
        else if (infoArray(b).is1D)  { return [...a, ...b];  }
        else if (infoArray(b).is2D)  { return [a].concat(b); }
    }
    else if (infoArray(a).is2D) {
        // Merge by adding rows if not specified
        if (dim == null) { dim = 0; }

        if (infoArray(b).is1D) {  
            if (dim == 0)  { return [...a].concat([b]); }
        }
        else if (infoArray(b).is2D) {
            if      (dim == 0)  { return [...a].concat(b);                                                      }
            else if (dim == 1)  { return new Array(a.length).fill().map((v,i) => [...a[i]].concat([...b[i]]) ); }
        }    
    }
}
function Cos(a) {
   let i, IA = infoArray(a);
   if      (isComplexNum(a)) { return a.Cos(); }
   else if (isNumber(a))     { return Math.cos(a); }
   else if (IA.is1D)         { return a.map((v,i) => Cos(v) ); }
   else if (IA.is2D)         { return a.map((v,i) => Cos(v) ); }
}
function Cosh(a) {
   let i, IA = infoArray(a);
   if      (isComplexNum(a)) { return a.Cosh(); }
   else if (isNumber(a))     { return Math.cosh(a); }
   else if (IA.is1D)         { return a.map((v,i) => Cosh(v) ); }
   else if (IA.is2D)         { return a.map((v,i) => Cosh(v) ); }
}
function Deg2Rad(a) {
   let i, IA = infoArray(a);
   if      (isComplexNum(a)) { return a.Deg2Rad(); }
   else if (isNumber(a))     { return a * Math.PI / 180; }
   else if (IA.is1D) { return a.map((v,i) => Deg2Rad(v) ); }
   else if (IA.is2D) { return a.map((v,i) => Deg2Rad(v) ); }
}
function Diag(a, k) {
   let N, IA = infoArray(a);
   if (IA.is1D) {
       // Returns a square matrix with the elements of Array-a[] array on the main diagonal.
       // k=0 represents the main diagonal, k>0 is above the main diagonal, and k<0 is below the main diagonal.
       if (k == null) { k = 0; };
       N = a.length;
       if   (k>=0) { return new Array(N+k).fill().map((vv,j) => new Array(N+k).fill().map((v,i) => {if(i==j+k) {return a[j];} else {return 0;}})); }
       else        { return new Array(N-k).fill().map((vv,j) => new Array(N-k).fill().map((v,i) => {if(i==j+k) {return a[i];} else {return 0;}})); }
   }
   else if (IA.is2D) {
       // Returns an array of the kth diagonal elements of Array-a[]
       // k=0 represents the main diagonal, k>0 is above the main diagonal, and k<0 is below the main diagonal.
       if (k == null) { k = 0; };
       N = a[0].length - Math.abs(k);
       if (N <= 0) { return 0; }
       else {
           if   (k >= 0) { return new Array(N).fill().map((v,i) => a[i][i+k]); }
           else          { return new Array(N).fill().map((v,i) => a[i-k][i]); }
       }
   }
}
function Diff(a, n) {
    // Differences and approximate derivatives
    // Calculates the nth difference of a-Array by applying the Diff(a) operator recursively n times 
    if (n==null) {n=1;}
    for (i=0; i<n; i++) { a = diff_(a); }
    return a;
    function diff_(a) {
        // Return if empty array
        if (a.length == 0) { return []; }
        let IA = infoArray(a);
        if      (isNumber(a)) { return []; }
        else if (isComplexNum(a)) { return []; }
        else if (IA.is1D) { if (a.length >= 2) { return Subtract(ExtractSubset(a, 1, a.length-1), ExtractSubset(a, 0, a.length-2)); } else { return []; } }
        else if (IA.is2D) { return a.map((v,i) => diff_(v)); }
    }
}
function Divide(a,b) {
   let IA = infoArray(a), IB = infoArray(b);
   if (isNumber(a)) {
       if      (isNumber(b))     { return a/b;                           }
       else if (isComplexNum(b)) { return new ComplexNum(a,0).Divide(b); }
       else if (IB.is1D)         { return b.map((v,i) => Divide(a, v) ); }
       else if (IB.is2D)         { return b.map((v,i) => Divide(a, v) ); }
   }
   else if (isComplexNum(a)) {
       if      (isNumber(b))     { return a.Divide(b); }
       else if (isComplexNum(b)) { return a.Divide(b); }
       else if (IB.is1D)         { return b.map((v,i) => Divide(a, v) ); }
       else if (IB.is2D)         { return b.map((v,i) => Divide(a, v) ); }
   }
   else if (IA.is1D) {
       if      (isNumber(b))     { return a.map((v,i) => Divide(v, b)    ); }
       else if (isComplexNum(b)) { return a.map((v,i) => Divide(v, b)    ); }
   }
   else if (IA.is2D) {
       if       (isNumber(b))    { return a.map((v,i) => Divide(a, b)); }
       else if (isComplexNum(b)) { return a.map((v,i) => Divide(a, b)); }
   }
}
function Dot(a, b) {
   // Dot product is a measure of how closely two vectors align, in terms of the directions they point.
   // if 2D-array, dot product is calculated for each row of a[,] and b[,] arrays
   let IA = infoArray(a), IB = infoArray(b);
   if      (isComplexNum(a) && isComplexNum(b)) { return a.Conj().Multiply(b); }
   else if (isComplexNum(a) && isNumber(b)    ) { return a.Conj().Multiply(b); }
   else if (isNumber(a)     && isComplexNum(b)) { return b.Multiply(a); }
   else if (isNumber(a)     && isNumber(b)    ) { return a*b;           }
   else if (IA.is1D && IB.is1D) { return Sum(a.map((v,i) => Dot(v, b[i]) )); }
   else if (IA.is2D && IB.is2D) { return a.map((v,i) => Dot(v, b[i]) ); }
}
function Eye(a, b) {
   let i, mm, temp=[];
   temp = Zeros(a, b);
   if (b == null) {
       temp[0] = 1;
       return temp;
   }
   else {
       mm   = Min([a, b]).val;
       for (i = 0; i <mm; i++) { temp[i][i] = 1; };
       return temp;
   }
}
function Exp(a) {
   let IA = infoArray(a);
   if      (isComplexNum(a)) { return a.Exp(); }
   else if (isNumber(a))     { return Math.exp(a); }
   else if (IA.is1D)         { return a.map((v,i) => Exp(v) ); }
   else if (IA.is2D)         { return a.map((v,i) => Exp(v) ); }
}
function ExtractSubset(a, a1, a2) {
    // A subset of data[] array is extracted between index numbers of a1 and a2.
    let IA = infoArray(a);
    if      (isNumber(a))     { return undefined; }
    else if (isComplexNum(a)) { return undefined; }
    else if (IA.is1D)         { return new Array(a2-a1+1).fill().map((v, i) => a[a1+i]); }
    else if (IA.is2D)         { return a.map((v,i) => ExtractSubset(v, a1, a2) ); }
}
function Flip(a) {
   let IA = infoArray(a);
   if      (isNumber(a)) { return a; }
   else if (isComplexNum(a)) { new ComplexNum(a.Re, a.Im); }
   else if (IA.is1D) { return [...a].reverse(); }
   else if (IA.is2D) { return a.map((v,i) => Flip(v)); }
}
function Find(a) {
    // Find indices of nonzero elements
    let IA = infoArray(a), res=[];
    if      (isNumber(a)) { if (a!=0) { return a; } else { return res; } }
    else if (isComplexNum(a)) { if ((a.Re!=0) && (a.Im!=0)) { return new ComplexNum(a.Re, a.Im);} else { return res; } }
    else if (IA.is1D) { a.map((v,i) => { if (v!=0) { res.push(i); } }); return res; }
    else if (IA.is2D) { a.map((v,i) => res.push(Find(v))); return res; }
}
function GetRange(a, rInd, cInd) {
    // rInd is an array containing the row index numbers
    // cInd is an array containing the column index numbers

    // Return empty array if a[] array is emty 
    if (IsArrayEmpty(a)) {return []; }

    // Check for dimesion of a[] array
    let IA = infoArray(a);
    
    // select all rows if not specified
    if (rInd == null) { rInd = LinStep(0, a.length-1, 1); }
 
    if      (IA.is1D) { return rInd.map((v, i) => a[v]) }
    else if (IA.is2D) {
        // select all columns if not specified
        if (cInd == null) { cInd = LinStep(0, a[0].length-1, 1); }
        return rInd.map((v, i) => GetRange(a[v], cInd) );
    }
 }
function GetColumn(a, ii) {
    let IA = infoArray(a), res=[];
    if      (IA.is1D) { return a[ii]; }
    else if (IA.is2D) { a.map((v,i) => res.push(GetColumn(v, ii))); return res; }
}
function IsContainNaN(a) {
   // Returns true if, at least, one of the elements of Array-a[] contains NaN
   // Runs for each row if 2D-Array-a[]
   let i=0, IA = infoArray(a);
   if      (isNumber(a))     { return Number.isNaN(a); }
   else if (isComplexNum(a)) { return a.IsNaN();       }
   else if (IA.is1D)         { while (i < a.length) { if (IsContainNaN(a[i])) { return true; }  i++; }; return false; }
   else if (IA.is2D)         { return a.map((v,i) => IsContainNaN(v) );  }
}
function IsArrayEmpty(a) {
    if (a.length == 0) { return true; } else { return false; }
}
function isComplexNum(a) {
   if (a.constructor.name.toUpperCase() === 'COMPLEXNUM') { return true; }
   else                                                   { return false; }
}
function IsEven(a) {
   // Returns TRUE if even number; otherwise, returns FALSE
   let IA = infoArray(a);
   if (isComplexNum(a)) { return false; }
   else if (isNumber(a)) { if (a%2 === 0) { return true; } else { return false; } }
   else if (IA.is1D) { return a.map((v,i) => IsEven(v)); }
   else if (IA.is2D) { return a.map((v,i) => IsEven(v)); }
}
function IsEqual(a, b) {
   let i, IA = infoArray(a), IB = infoArray(b);
   if      (isComplexNum(a) && isComplexNum(b)) { return a.IsEqual(b); }
   else if (isComplexNum(a) && isNumber(b)    ) { return a.IsEqual(b); }
   else if (isNumber(a)     && isComplexNum(b)) { return b.IsEqual(a); }
   else if (isNumber(a)     && isNumber(b)    ) { return a === b;      }
   else if (IA.is1D) {
       if      (isComplexNum(b) || isNumber(b)) { return a.map((v, i) => IsEqual(v, b   ) ); }
       else if (IB.is1D)                        { return a.map((v, i) => IsEqual(v, b[i]) ); }
   }
   else if (IA.is2D) {
       if      (isComplexNum(b) || isNumber(b)) { return a.map((v,i) => IsEqual(v, b   ) ); }
       else if (IB.is2D)                        { return a.map((v,i) => IsEqual(v, b[i]) ); }
   }
}
function IsFinite(a) {
    // returns a logical array of the same size as the input, where true indicates 
    // a finite value and false indicates a non-finite value (i.e., Inf, -Inf, or NaN)

    if (IsArrayEmpty(a)) { return 0; }

    let IA = infoArray(a);
    if (isComplexNum(a))  { return Number.isFinite(a.Re) && Number.isFinite(a.Im) ? true: false; }
    else if (isNumber(a)) { return Number.isFinite(a) ? true: false;  }
    else if (IA.is1D)     { return a.map((v,i) => IsFinite(v)); }
    else if (IA.is2D)     { return a.map((v,i) => IsFinite(v)); }
}
function IsNaN(a) {
   let IA = infoArray(a);
   if      (isComplexNum(a)) { return a.IsNaN(); }
   else if (isNumber(a))     { return Number.isNaN(a); }
   else if (IA.is1D)         { return a.map((v,i) => IsNaN(v) ); }
   else if (IA.is2D)         { return a.map((v,i) => IsNaN(v) ); }
}
function isNumber(a) {
   if (a.constructor.name.toUpperCase() === 'NUMBER') { return true; }
   else                                               { return false; }
}
function IsPosDef(a) {
    // Return true, if 2D=Array-a[] is positive definite
    if (Chol(a).isPosDef) { return true; } else { return false; }
}
function IsSymmetric(a) {
   // Returns true if Array-a[] is symmetric; otherwise, returns false.
   var i, j, aa, bb, IA = infoArray(a);
   if      (isNumber(a) || isComplexNum(a)) { return true;  }
   else if (IA.is1D)                        { return false; }
   else if (IA.is2D) {
       for (i=0; i<a.length; i++) {
           for (j=i+1; j<a[0].length; j++) {
               aa = a[i][j];
               bb = a[j][i];
               if      (isComplexNum(aa) && isComplexNum(bb)) { if (!aa.IsEqual(bb)) { return false; } }
               else if (isComplexNum(aa) && isNumber(bb)    ) { if (!aa.IsEqual(bb)) { return false; } }
               else if (isNumber(aa)     && isComplexNum(bb)) { if (!bb.IsEqual(aa)) { return false; } }
               else if (isNumber(aa)     && isNumber(bb)    ) { if (a[i][j] != a[j][i])  { return false; } }
           }
       }
       return true;
   }
   else { return false; }
}
function Imag(a) {
   let i, IA = infoArray(a);
   if      (isComplexNum(a)) { return a.Imag(); }
   else if (isNumber(a))     { return 0;        }
   else if (IA.is1D)         { return a.map((v,i) => Imag(v) ); }
   else if (IA.is2D)         { return a.map((v,i) => Imag(v) ); }
}
function Inv(a) {
    // Returns the inverse of square 2D-Array-a[].

   let n, m, IA=[], Re, i, j, k, x=[], y=[], sum, IAA = infoArray(a);

   if      (IsEqual(1/a, Number.POSITIVE_INFINITY)) { return Number.POSITIVE_INFINITY; }
   else if (IsEqual(1/a, Number.NEGATIVE_INFINITY)) { return Number.NEGATIVE_INFINITY; }
   else if (isComplexNum(a))                        { return a.Pow(-1); }
   else if (isNumber(a))                            { return 1/a; }

   // Return if not 2D-Array-a[]
   if (!IAA.is2D) { throw new Error("Matrix must be square."); }

   n  = a.length;
   m  = a[0].length
   IA = Zeros(n, n);

   // Return if not square
   if (n!=m) { throw new Error("Matrix must be square."); }

   // Calculate LU decomposition such that P*a = L*U
   Re = LU(a);  // L = Re.L;  U = Re.U;  P = Re.P;

   for (j=0; j<n; j++) {
       // Find solution of L*y = B
       y = new Array(n);
       x = new Array(n);

       for (i=0; i<n; i++) {
           sum = 0;
           for (k=0; k<i; k++) {
               sum = Add(sum, Mult( Re.L[i][k], y[k]));
           }
           y[i] = Subtract(Re.P[i][j], sum);
       }

       // find solution of U*x = y
       for (i=n-1; i>=0; i--) {
           sum = 0;
           for (k=i+1; k<n; k++) {
               sum = Add(sum, Mult(Re.U[i][k], x[k]));
           }
           x[i] = Mult(Divide(1, Re.U[i][i]), Subtract(y[i], sum));
       }

       // x[] array is the column of the inverse IA[,]
       for (i=0; i<n; i++) {
           IA[i][j] = x[i];
       }
   }
   return IA;
}
function LinSpace(x1, x2, N) {
   // Returns an array of N-points evenly spaced points between x1 and x2.
   let i, y, step;
   if (isNumber(x1) && isNumber(x2)) {
       if (N <= 0) {
           // Return empty array
           return new Array(0).fill();
       }
       else if (N == 1) {
           // Return x2
           return [x2];
       }
       else if (N > 1) {
           if (x1 === x2) { return new Array(N).fill(x1) }
           y    = new Array(N).fill();
           step = (x2 - x1) / (N - 1);
           y[0] = x1;
           for (i = 1; i < N-1; i++) { y[i] = y[i - 1] + step; }
           y[N - 1] = x2;
           return y;
       }
   }
   else if (isNumber(x1)     && isComplexNum(x2)) { return Num2Complex(LinSpace(x1,    x2.Re, N),   LinSpace(0,     x2.Im, N)); }
   else if (isComplexNum(x1) && isNumber(x2))     { return Num2Complex(LinSpace(x1.Re, x2,    N),   LinSpace(x1.Im, 0,     N)); }
   else if (isComplexNum(x1) && isComplexNum(x2)) { return Num2Complex(LinSpace(x1.Re, x2.Re, N),   LinSpace(x1.Im, x2.Im, N)); }
}
function LinStep(x1, x2, step) {
   // Return an array with evenly spaced points between x1 and x2. Array values may exceed x2.
   // Number of points depends on the step value.
   let i, N, y, alfa, Dx, Dy, N1, N2;

   // Always work with absolute value of step
   step = Abs(step);

   if (isNumber(x1) && isNumber(x2)) {
       if      (x2 === x1) { return [x1]; }
       else if (x1 > x2)   { N = Math.ceil( (x1 - x2) / step ) + 1;  step *= -1; }
       else if (x2 > x1)   { N = Math.ceil( (x2 - x1) / step ) + 1;              }

       y = new Array(N).fill();
       y[0] = x1;
       for (i = 1; i < Abs(N); i++) { y[i] = y[i-1] + step }
       return y;
   }
   else if (isNumber(x1) && isComplexNum(x2)) {
       Dx   = x2.Re - x1;
       Dy   = x2.Im - 0;
       alfa = Math.atan( Dy / Dx );
       N1   = step * Math.cos(alfa);
       N2   = step * Math.sin(alfa);
       if      ((Dx !== 0) && (Dy !== 0)) { return Num2Complex(LinStep(x1, x2.Re, N1),   LinStep(0, x2.Im, N2)); }
       else if ((Dx === 0) && (Dy !== 0)) { return Num2Complex(x1,                       LinStep(0, x2.Im, N2)); }
       else if ((Dx !== 0) && (Dy === 0)) { return Num2Complex(LinStep(x1, x2.Re, N1),   x2.Im                ); }
       else if ((Dx === 0) && (Dy === 0)) { return Num2Complex(x1,                       x2.Im                ); }
   }
   else if (isComplexNum(x1) && isNumber(x2)) {
       Dx   = x2 - x1.Re;
       Dy   = 0 - x1.Im;
       alfa = Math.atan( Dy / Dx );
       N1   = step * Math.cos(alfa);
       N2   = step * Math.sin(alfa);
       if      ((Dx !== 0) && (Dy !== 0)) { return Num2Complex(LinStep(x1.Re, x2, N1),   LinStep(x1.Im, 0, N2)); }
       else if ((Dx === 0) && (Dy !== 0)) { return Num2Complex(x1.Re,                    LinStep(x1.Im, 0, N2)); }
       else if ((Dx !== 0) && (Dy === 0)) { return Num2Complex(LinStep(x1.Re, x2, N1),   0                    ); }
       else if ((Dx === 0) && (Dy === 0)) { return Num2Complex(x1.Re,                    0                    ); }
   }
   else if (isComplexNum(x1) && isComplexNum(x2)) {
       Dx   = x2.Re - x1.Re;
       Dy   = x2.Im - x1.Im;
       alfa = Math.atan( Dy / Dx );
       N1   = step * Math.cos(alfa);
       N2   = step * Math.sin(alfa);
       if      ((Dx !== 0) && (Dy !== 0)) { return Num2Complex(LinStep(x1.Re, x2.Re, N1),   LinStep(x1.Im, x2.Im, N2)); }
       else if ((Dx === 0) && (Dy !== 0)) { return Num2Complex(x1.Re,                       LinStep(x1.Im, x2.Im, N2)); }
       else if ((Dx !== 0) && (Dy === 0)) { return Num2Complex(LinStep(x1.Re, x2.Re, N1),   x2.Im                    ); }
       else if ((Dx === 0) && (Dy === 0)) { return Num2Complex(x1.Re,                       x2.Im                    ); }
   }
}
function Log(a) {
   let i, IA = infoArray(a);
   if      (isComplexNum(a)) { return a.Log(); }
   else if (isNumber(a))     { return Math.log(a); }
   else if (IA.is1D)         { return a.map((v,i) => Log(v) ); }
   else if (IA.is2D)         { return a.map((v,i) => Log(v) ); }
}
function Log10(a) {
   let i, IA = infoArray(a);
   if      (isComplexNum(a)) { return a.Log10(); }
   else if (isNumber(a))     { return Math.log10(a); }
   else if (IA.is1D)         { return a.map((v,i) => Log10(v) ); }
   else if (IA.is2D)         { return a.map((v,i) => Log10(v) ); }
}
function LU(a, Opt) {
   // Returns the L[], U[] and P[] arrays such that P*a = L*U
   // If Opt is true, the algorithm will use pivots (row exchange)
   // If Opt is set to false, no pivoting will be used unless the leading number is equal to or very close to zero.
   // In such case, pivoting must (will) be used to eliminate the numerical instability.
   let i, j, k, n, m, t, P=[], L=[], U=[], UU=[], VM, IA=infoArray(a);

   if (isNumber(a) || isComplexNum(a) || (IA.is1D)) { return {'L': 1, 'U': a, 'P': 1 } }

   if (Opt == null) { Opt = true; }

   n = a.length;
   m = a[0].length;
   t = Math.min(n, m);
   L = Zeros(n, t);
   U = Copy(a);
   P = Eye(n, n);

   //This will calculate the L[] and U[] arrays with optional row exchange
   for (i=0; i<t; i++) {
       // Row-exchange
       if ( IsEqual(U[i][i], 0) || Opt)  {
           // Determine the maximum value and the Index number of each column of the U[,] array
           VM = Max(Abs(GetRange(U, LinStep(i, n-1, 1), [i]).flat())).Indx;
           SwapRow(U, [VM + i], [i])
           SwapRow(P, [VM + i], [i])
           SwapRow(L, [VM + i], [i])
       }
       for (j=i+1; j<n; j++) {
           L[j][i] = Divide( U[j][i], U[i][i] );
           for (k=i; k<m; k++) {
               //for each element in the row
               U[j][k] = Subtract( U[j][k], Mult( L[j][i], U[i][k] ) );
           }
       }
       L[i][i] = 1;
   }

   if (t < n) {
       // Resize the U[,] matrix
       UU = Zeros(t, m);
       for (i=0; i<t; i++) { for (j=0; j<m; j++) { UU[i][j] = U[i][j]; }}
       return {'L': L, 'U': UU, 'P': P }
   }
   else {
       return {'L': L, 'U': U,  'P': P }
   }
}
function Max(a) {
   let i, maxVal, maxInd, temp, IA = infoArray(a);
   if      (isNumber(a))     { return { 'val': a, 'Indx': undefined }; }
   else if (isComplexNum(a)) { return { 'val': new ComplexNum(a.Re, a.Im), 'Indx': undefined }; }
   else if (IA.is1D) {
       maxVal = Number.NEGATIVE_INFINITY;
       maxInd = undefined;
       a.map((v,i) => { if (isComplexNum(v)) { temp = Abs(v); } else {temp=v; }; if (temp > maxVal) { maxVal = temp; maxInd = i; } } )
       temp = a[maxInd];
       if      (isComplexNum(temp)) { return { 'val': new ComplexNum(temp.Re, temp.Im), 'Indx': maxInd }   }
       else if (isNumber(temp))     { return { 'val': temp,                             'Indx': maxInd } }
   }
   else if (IA.is2D) {
       maxVal = [];
       maxInd = [];
       for (i = 0; i < a.length; i++) {
            temp = Max(a[i]);
            maxVal.push(temp.val);
            maxInd.push(temp.Indx);
       }
       return { 'val': maxVal, 'Indx': maxInd }
   }
}
function Min(a) {
   let i, minVal, minInd, temp, IA = infoArray(a);
   if      (isNumber(a))     { return { 'val': a, 'Indx': undefined }; }
   else if (isComplexNum(a)) { return { 'val': new ComplexNum(a.Re, a.Im), 'Indx': undefined }; }
   else if (IA.is1D) {
       minVal = Number.POSITIVE_INFINITY;
       minInd = undefined;
       a.map((v,i) => { if (isComplexNum(v)) { temp = Abs(v); } else {temp=v; }; if (temp < minVal) { minVal = temp; minInd = i; } } )
       temp = a[minInd];
       if      (isComplexNum(temp)) { return { 'val': new ComplexNum(temp.Re, temp.Im), 'Indx': minInd }   }
       else if (isNumber(temp))     { return { 'val': temp,                             'Indx': minInd } }
   }
   else if (IA.is2D) {
       minVal = [];
       minInd = [];
       for (i = 0; i < a.length; i++) {
           temp = Min(a[i]);
           minVal.push(temp.val);
           minInd.push(temp.Indx);
       }
       return { 'val': minVal, 'Indx': minInd }
   }
}
function Mult(a, b) {

    if (isNumber(a)) {
        if      (IsArrayEmpty(b))    { return [];                           }
        else if (isNumber(b))        { return a * b;                        }
        else if (isComplexNum(b))    { return b.Multiply(a);                }
        else if (infoArray(b).is1D)  { return b.map((v,i) => Mult(a, v));   }
        else if (infoArray(b).is2D)  { return b.map((v,i) => Mult(a, v));   }
    }
    else if (isComplexNum(a)) {
        if      (IsArrayEmpty(b))    { return [];                           }
        else if (isNumber(b))        { return a.Multiply(b);                }
        else if (isComplexNum(b))    { return a.Multiply(b);                }
        else if (infoArray(b).is1D)  { return b.map((v,i) => Mult(a, v));   }
        else if (infoArray(b).is2D)  { return b.map((v,i) => Mult(a, v));   }
    }
    else if (IsArrayEmpty(a)) {
        if      (IsArrayEmpty(b))    { return [];    }
        else if (isNumber(b))        { return [];    }
        else if (isComplexNum(b))    { return [];    }
    }
    else if (infoArray(a).is1D) { 
        if      (isNumber(b))        { return a.map((v,i) => v*b);                                                }
        else if (isComplexNum(b))    { return a.map((v,i) => b.Multiply(v));                                      }
        else if (infoArray(b).is2D)  { return b[0].map((v,i) => Sum(Mult_El(a, GetRange(b, null, [i]).flat())));  }
    }
    else if (infoArray(a).is2D) {
        if      (isNumber(b))        { return a.map((v,i) => Mult(v, b)); }
        else if (isComplexNum(b))    { return a.map((v,i) => Mult(v, b)); }
        else if (infoArray(b).is2D)  { return a.map((v,i) => Mult(v, b)); }
    }
    else {
        throw new Error("ERROR IN MATRIX MULTIPLICATION - MATRIX DIMENSIONS DO NOT MATCH");
    }
}
function Mult_old(a, b) {

    // Return empty array if a[]-array is empty 
    if (IsArrayEmpty(a)) { return []; }

    let IA = infoArray(a), IB = infoArray(b);
    if      (isNumber(a)     && isNumber(b))     { return a * b;         }
    else if (isNumber(a)     && isComplexNum(b)) { return b.Multiply(a); }
    else if (isComplexNum(a) && isNumber(b))     { return a.Multiply(b); }
    else if (isComplexNum(a) && isComplexNum(b)) { return a.Multiply(b); }
    else if (isNumber(a) || isComplexNum(a)) {
        if      (IB.is1D) { return b.map((v,i) => Mult(a, v) ); }
        else if (IB.is2D) { return b.map((v,i) => Mult(a, v) ); }
    }
    else if (isNumber(b) || isComplexNum(b)) {
        if      (IA.is1D) { return a.map((v,i) => Mult(v, b) ); }
        else if (IA.is2D) { return a.map((v,i) => Mult(v, b) ); }
    }
    else if (IA.is1D && IB.is2D) { return b[0].map((v,i) => Sum(Mult_El(a, GetRange(b, null, [i]).flat())) ); }      // a=[1x4]  b=[4x2]   OR   a=[1x4]  b=[4x1]
    else if (IA.is2D && IB.is2D) { return a.map((v,i) => Mult(v, b) ); }                                             // a=[2x3]  b=[3x4]   OR   a=[2x1]  b=[1x4] 
    else if (IA.is2D && IB.is1D) { return a.map((v,i) => Sum(Mult_El(v, b))); }                                      // a=[2x3]  b=[3x1]
}
function Mult_El(a, b) {
   // Element-wise multiplication ONLY
   let IA = infoArray(a), IB = infoArray(b);
   if      (isNumber(a)     && isNumber(b))     { return a * b;         }
   else if (isNumber(a)     && isComplexNum(b)) { return b.Multiply(a); }
   else if (isComplexNum(a) && isNumber(b))     { return a.Multiply(b); }
   else if (isComplexNum(a) && isComplexNum(b)) { return a.Multiply(b); }
   else if (IA.is1D && IB.is1D)                 { return a.map((v,i) => Mult_El(a[i], b[i]) ); }
   else if (IA.is2D && IB.is2D)                 { return a.map((v,i) => Mult_El(v,    b[i]) ); }
}
function Mean(a) {
   let temp, IA = infoArray(a);
   if      (isNumber(a))     { return a;                          }
   else if (isComplexNum(a)) { return new ComplexNum(a.Re, a.Im); }
   else if (IA.is1D)         {
       temp = Sum(a);
       if      (isNumber(temp))     { return temp/a.length;         }
       else if (isComplexNum(temp)) { return temp.Divide(a.length); }
   }
   else if (IA.is2D) { return a.map((v,i) => Mean(v) ); }
}
function Median(a) {
    let temp, N, IA = infoArray(a);
    if      (isNumber(a))     { return a;                          }
    else if (isComplexNum(a)) { return new ComplexNum(a.Re, a.Im); }
    else if (IA.is1D)         {
        temp = Sort(a, "ASC").data;
        N    = temp.length;
        if (IsEven(N)) { return (temp[N/2]+temp[N/2-1])/2; } else { return temp[(N-1)/2]; }
    }
    else if (IA.is2D) { return a.map((v,i) => Meadian(v) ); }
}
function NextPow2(a) {
   // Returns the next unsigned integer power of 2
   let i, power=1, IA = infoArray(a);
   if      (isNumber(a))     { while (power < a) { power *= 2; }; return power; }
   else if (isComplexNum(a)) { return NextPow2( Abs(a) );                       }
   else if (IA.is1D)         { return a.map((v,i) => NextPow2(v) );             }
   else if (IA.is2D)         { return a.map((v,i) => NextPow2(v) );             }
}
function Norm(a, p) {
    // Returns the generalized vector p-norm.
    // If p is not given, Eucliden norm is returned. (The p=2 norm is also called the 2-norm, vector magnitude, or Euclidean length)
    if (IsArrayEmpty(a)) { return 0; }

    let IA = infoArray(a);
    if (p == null) { p = 2; }

    if      (isNumber(a))     { return a;      }
    else if (isComplexNum(a)) { return Abs(a); }
    else if (IA.is1D) {
        if      (p == Number.POSITIVE_INFINITY)  { return Max(Abs(a)).val;                    }
        else if (p == Number.NEGATIVE_INFINITY)  { return Min(Abs(a)).val;                    }
        else if (p == 0)                         { return Number.POSITIVE_INFINITY;           }
        else if (p == 1)                         { return Sum( Abs(a) );                      }
        else if ((p > 0) && (isNumber(p)))       { return Math.pow(Sum(Pow(Abs(a), p)), 1/p); }
        else                                     { return undefined; }

    }
    else if (IA.is2D) {
        if      (p == Number.POSITIVE_INFINITY)  { return Max(Sum(Abs(Transpose(a)))).val;     }
        else if (p == Number.NEGATIVE_INFINITY)  { return undefined;                           }
        else if (p == 1)                         { return Max(a.map((v,i) => Norm(v, p))).val; }
        else if (p == 2)                         { return Sqrt(Max(Matrix_Eig(Mult(Transpose(A), A))[0]).val); }  //  sqrt(max(eig(a' * a))), but SVD function is missing       else if ((p > 0) && (isNumber(p)))       { return undefined                            }
        else if (p === 'fro')                    { return Math.sqrt( Sum(Sum(Pow(Abs(a),2))) ) }
    }
}
function Num2Complex(a, b) {
   // a=> Real part
   // b=> Imaginary part
   let i, IA = infoArray(a), IB = infoArray(b);

   if (isNumber(a)) {
       if      (b == null)   { return new ComplexNum(a, 0); }
       else if (isNumber(b)) { return new ComplexNum(a, b); }
       else if (IB.is1D)     { return b.map((v,i) => new ComplexNum(a, v) ); }
   }
   else if (IA.is1D) {
       if      (b == null)   { return a.map((v,i) => new ComplexNum(v, 0) );    }
       else if (isNumber(b)) { return a.map((v,i) => new ComplexNum(v, b) );    }
       else if (IB.is1D)     { return a.map((v,i) => new ComplexNum(v, b[i]) ); }
   }
   else if (IA.is2D) {
       if      (b == null)   { return a.map((v,i) => Num2Complex(v      ) ); }
       else if (isNumber(b)) { return a.map((v,i) => Num2Complex(v, b   ) ); }
       else if (IB.is2D)     { return a.map((v,i) => Num2Complex(v, b[i]) ); }
   }
}
function Ones(a, b) {
   if      (b == null) { return new Array(a).fill(1); }
   else if (b != null) { return new Array(a).fill().map(() => Ones(b)); }
}
function Percentile(a, pr) {
    // Returns percentiles of elements in a[] array for the percentages pr in the interval [0, 100]
    // Return percentile for each array if a[] is 2D array
    // Percentile treats NaN valies as missing values and removes them from a[]
    let temp, N, R, FR, IA = infoArray(a);
    if      (isNumber(a))     { return a;                          }
    else if (isComplexNum(a)) { return new ComplexNum(a.Re, a.Im); }
    else if (IA.is1D)         {
        // Remove NaN and sort a[] array in ascending order
        temp = Sort(RemoveNaN(a), "ASC").data;
        N = temp.length;
        if (infoArray(pr).is1D) {
            let result = [];
            pr.map((v,i) => {
                if (N == 1) { result.push(a[0]); }
                else {
                    R = (v/100) * (N - 1);
                    FR = R - Math.trunc(R);
                    R = Math.trunc(R);
                    result.push(FR * (temp[R+1] - temp[R]) + temp[R]);
                }
            });
            return result;
        }
        else {
            if (N == 1) { return a[0]; }
            else {
                R = (pr/100) * (N - 1);
                FR = R - Math.trunc(R);
                R = Math.trunc(R);
                return FR * (temp[R+1] - temp[R]) + temp[R];
            }
        }
    }
    else if (IA.is2D) { return a.map((v,i) => Percentile(v, pr) ); }
}
function Pow(a, p) {
    let T, IA = infoArray(a);
    if      (isNumber(a))     { return Math.pow(a, p); }
    else if (isComplexNum(a)) { return a.Pow(p);       }
    else if (IA.is1D)         { return a.map((v,i) => Pow(v,p) ); }
    else if (IA.is2D)         { return a.map((v,i) => Pow(v,p) ); }
}
function Prod(a) {
    let IA = infoArray(a);
    if      (isNumber(a))     { return a; }
    else if (isComplexNum(a)) { return a; }
    else if (IA.is1D)         { return a.reduce((a, num) => Mult(a, num), 1);  }
    else if (IA.is2D)         { return a.map((v,i) => Prod(v,p) ); }

}
function Print(a, N, M) {
   let i, m,n, Str=[], IA = infoArray(a);

   if (N==null) { N=3; }
   if (M==null) { M=4; }

   if (IA.is1D) {
       n = a.length;
       m = Math.min(n, N);
       for (i=0; i<m; i++) {
           R = Real(a[i]);
           I = Imag(a[i]);
           if (Sign(I)>=0) { ss = '+'; } else { ss = '-'; }
           Str += R.toFixed(M) + ss + Abs(I).toFixed(M)+'i  ';
       }
       Str = Str.trim();
       if (n>=2*N) {
           m = Math.max(N, n-N);
           Str += '  ...........  ';

           for (i=m; i<n; i++) {
               R = Real(a[i]);
               I = Imag(a[i]);
               if (Sign(I)>=0) { ss = '+'; } else { ss = '-'; }
               Str += R.toFixed(M) + ss + Abs(I).toFixed(M)+'i ';
           }
           Str = Str.trim();
       }
   }
   if (IA.is2D) {
       a.map((v,i) => Str.push(Print(v, N, M))  )
   }
}
function QR(a, Option) {
   // This function uses Householder reflections to decompose 2D-Array-a[] into a product of a = Q * R
   //
   // A Householder reflection (or Householder transformation) is a transformation that takes a vector and reflects it
   // about some plane or hyperplane to create a mirror copy.
   // We can use this operation to calculate the QR factorization of an m-by-n 2D-Array-a[] with with m ≥ n.
   //
   // QR decomposition/factorization is often used  to solve the linear least squares problem and is also the basis for
   // a particular eigenvalue algorithm, the QR algorithm.
   //
   // A[,] is a rectangular matrix of m-by-n size
   // Q[,] is an orthonormal, unitary matrix of m-by-m size
   // R[,] is an upper triangular matrix of m-by-n size
   //
   // If Option = 1, then it produces an economy-size decomposition.
   // The size of the outputs depends on the size of m-by-n matrix A[,]
   //
   // This function is cross-checked against Matlab [Q,R]=qr(A) built-in function 2023-06-08

   // Check default values of input arguments
   if (Option == null) { Option = 0; };

   let m, n, R=[], Q=[], i, j, k, tt, temp, Vnorm, s, u1, u, w, tau, aa, IndR, IndC;
   let IA = infoArray(a);

   if      (isNumber(a))     { return {'Q': 1, 'R': a}; }
   else if (isComplexNum(a)) { return {'Q': Divide(a, Abs(a)), 'R': Abs(a)}; }
   else if (IA.is1D)         { let Re = QR([a]); return {'Q': Re.Q[0][0], 'R': Re.R[0]}; }
   else if (IA.is2D) {
       m = a.length;
       n = a[0].length;
       R = Copy(a);
       Q = Eye(m, m);

       if (m >= n) { tt = n; } else { tt = m; }

       for (i=0; i < tt; i++) {

           temp  = GetRange(R, LinStep(i, m-1, 1), [i]).flat();
           Vnorm = Norm(temp);
           s     = Mult(-1, Sign(temp[0]));
           u1    = Subtract(temp[0], Mult(s, Vnorm));
           u     = Copy(temp); u[0] = u1;
           w     = Divide(u, u1);
           tau   = Mult(Mult(-1,Conj(s)), Divide(u1, Vnorm));

           // Update R[,] matrix, one column at a time
           for (k=0; k<n; k++) {
               aa = Mult(tau, Dot(w, GetRange(R, LinStep(i, m-1, 1), [k]).flat()));
               for (j=i; j<m; j++) { R[j][k] = Subtract(R[j][k], Mult(aa, w[j-i])); }
           }

           // Update Q[,] matrix: one row at a time
           for (k = 0; k<m; k++) {
               aa = Mult(tau, Sum(Mult_El(w, GetRange(Q, [k], LinStep(i, m-1, 1)).flat())));
               for (j=i; j<m; j++) { Q[k][j] = Subtract(Q[k][j], Mult(aa, Conj(w[j-i]))) }
           }
       }
       if ((Option == 1) && (m > n)) {
           // Return economy-size decomposition such that
           // [Q] m-by-n   and  [R] n-by-n
           IndR = LinStep(0, m-1, 1);
           IndC = LinStep(0, n-1, 1);
           return {'Q': GetRange(Q, IndR, IndC), 'R': GetRange(R, IndC, IndC) }
       } else { return {'Q': Q, 'R': R}; }
   }

}
function Rand(a, b) {
   let IA = infoArray(a), IB = infoArray(b);
   if (a==null) {
       if      (b==null)            { return Math.random(); }
       else if (isNumber(b))        { throw new Error("Inconsistent Arguments."         ); }
       else if (isComplexNum(b))    { throw new Error("Arguments must be real numbers." ); }
       else if (IB.is1D || IB.is2D) { throw new Error("Arguments cannot be array."      ); }
   }
   else if (isNumber(a)) {
       if      (b==null)            { return new Array(a).fill().map(() => Math.random()); }
       else if (isNumber(b))        { return new Array(a).fill().map(() => Rand(b)      ); }
       else if (isComplexNum(b))    { throw new Error("Arguments must be real numbers." ); }
       else if (IB.is1D || IB.is2D) { throw new Error("Arguments cannot be array."      ); }
   }
   else if (isComplexNum(a))        { throw new Error("Arguments must be real numbers." ); }
   else if (IA.is1D || IA.is2D)     { throw new Error("Arguments cannot be array."      ); }
}
function Rand_Complex(a, b) {
   let IA = infoArray(a), IB = infoArray(b);
   if (a==null) {
       if      (b==null)            { return new ComplexNum(Math.random(), Math.random()); }
       else if (isNumber(b))        { throw new Error("Arguments ar inconsistent."); }
       else if (isComplexNum(b))    { throw new Error("Arguments must be real numbers."); }
       else if (IB.is1D || IB.is2D) { throw new Error("Arguments cannot be array."); }
   }
   else if (isNumber(a)) {
       if      (b==null)            { return new Array(a).fill().map(() => new ComplexNum(Math.random(), Math.random())); }
       else if (isNumber(b))        { return new Array(a).fill().map(() => Rand_Complex(b)                             ); }
       else if (isComplexNum(b))    { throw new Error("Arguments must be real numbers."); }
       else if (IB.is1D || IB.is2D) { throw new Error("Arguments cannot be array."); }
   }
   else if (isComplexNum(a))        { throw new Error("Arguments must be real numbers."); }
   else if (IA.is1D || IA.is2D)     { throw new Error("Arguments cannot be array."); }
}
function Randn(a, b) {
   // Generates a N-length array of random numbers that are normally distributed
   // Normal distribution =>   mu = 0;  std = 1;
   d = 2 * Math.PI;
   c = Math.sqrt(2 * Math.PI);
   if (b == null) {
       return new Array(a).fill().map(() => {
           U1 = Math.random();
           U2 = Math.random();
           if (U1 != 0) { x = Math.sqrt(-2 * Math.log(U1)) * Math.cos(d * U2); return Math.exp(-x * x / 2) / c; } else { return 0; }
       });
   }
   else if (b != null) { return new Array(a).fill().map(() => Randn(b) ); }
}
function Real(a) {
   let IA = infoArray(a);
   if      (isComplexNum(a)) { return a.Real(); }
   else if (isNumber(a))     { return a;        }
   else if (IA.is1D)         { return a.map((v,i) => Real(v) ); }
   else if (IA.is2D)         { return a.map((v,i) => Real(v) ); }
}
function RemoveNaN(a) {
    // Removes elements of NaN
    let result=[], IA = infoArray(a);
    if (IA.is1D) { 
        IsNaN(a).map((v,i) => {if(!v){result.push(i);}}); 
        return GetRange(a, result);
    }
    else if (IA.is2D) { return a.map((v,i) =>  RemoveNaN(v) ); }
}
function ReplaceRow(a, Ind, newVal) {
    // Replaces rows of 2D-a[]-Array with the corresponding rows of newVal array
    // Ind is contains the row index-numbers of a[] array - could be an array of row numbers or a single row number 
    // newVal is the new array to be replaced 
    let IA = infoArray(Ind);
    if (isNumber(Ind)) { a[Ind] = newVal; }
    else if (IA.is1D)  { Ind.map((v,i) => ReplaceRow(a, v, newVal[i])); }
    return a;
}
function ReplaceCol(a, Ind, newVal) {
    // Replaces columns of 2D-a[]-Array with the corresponding columns of newVal array
    // Ind is an array containing the column index numbers
    // newVal is an array containing the new vales
    let i, IA = infoArray(Ind);
    if (isNumber(Ind)) { for (i=0; i<a.length; i++) { a[i][Ind] = newVal[i]; }  }
    else if (IA.is1D)  { Ind.map((v, ii) => ReplaceCol(a, v, newVal[ii])); }
    return a;
}
function Rms(a) {
   let i, IA = infoArray(a);
   if      (isNumber(a))     { return a;       }
   else if (isComplexNum(a)) { return Abs(a);  }
   else if (IA.is1D)         { return Math.sqrt(Mean(Pow(Abs(a), 2))); }
   else if (IA.is2D)         { return a.map((v,i) => Rms(v));  }
}
function Sign(a) {
   let i, IA = infoArray(a);
   if      (isComplexNum(a)) { return a.Sign(); }
   else if (isNumber(a))     { return Math.sign(a); }
   else if (IA.is1D)         { return a.map((v,i) => Sign(v) ); }
   else if (IA.is2D)         { return a.map((v,i) => Sign(v) ); }
}
function Sin(a) {
   let i, IA = infoArray(a);
   if      (isComplexNum(a)) { return a.Sin(); }
   else if (isNumber(a))     { return Math.sin(a); }
   else if (IA.is1D)         { return a.map((v,i) => Sin(v) ); }
   else if (IA.is2D)         { return a.map((v,i) => Sin(v) ); }
}
function Sinh(a) {
   let i, IA = infoArray(a);
   if      (isComplexNum(a)) { return a.Sinh(); }
   else if (isNumber(a))     { return Math.sinh(a); }
   else if (IA.is1D)         { return a.map((v,i) => Sinh(v) ); }
   else if (IA.is2D)         { return a.map((v,i) => Sinh(v) ); }
}
function Std(a) {
    let IA = infoArray(a);
   if      (isComplexNum(a)) { return 0; }
   else if (isNumber(a))     { return 0; }
   else if (IA.is1D)         { return Math.sqrt(Var(a));       }
   else if (IA.is2D)         { return a.map((v,i) => Std(v) ); }
}
function Subtract(a, b) {
   let i, IA = infoArray(a), IB = infoArray(b);
   if      (isNumber(a)) {
       if      (isNumber(b))     { return a - b;                          }
       else if (isComplexNum(b)) { return new ComplexNum(a-b.Re, -b.Im);  }
       else if (IB.is1D)         { return b.map((v,i) => Subtract(a, v)); }
       else if (IB.is2D)         { return b.map((v,i) => Subtract(a, v)); }
   }
   else if (isComplexNum(a)) {
       if      (isNumber(b))     { return a.Subtract(b); }
       else if (isComplexNum(b)) { return a.Subtract(b); }
       else if (IB.is1D)         { return b.map((v,i) => Subtract(a, v)); }
       else if (IB.is2D)         { return b.map((v,i) => Subtract(a, v)); }
   }
   else if (IA.is1D) {
       if      (isNumber(b))     { return a.map((v,i) => Subtract(v, b   )); }
       else if (isComplexNum(b)) { return a.map((v,i) => Subtract(v, b   )); }
       else if (IB.is1D)         { return a.map((v,i) => Subtract(v, b[i])); }
       else if (IB.is2D)         { return b.map((v,i) => Subtract(a, v   )); }
   }
   else if (IA.is2D) {
       if      (isNumber(b))     { return a.map((v,i) => Subtract(v, b   )); }
       else if (isComplexNum(b)) { return a.map((v,i) => Subtract(v, b   )); }
       else if (IB.is1D)         { return a.map((v,i) => Subtract(v, b   )); }
       else if (IB.is2D)         { return a.map((v,i) => Subtract(v, b[i])); }
   }
}
function Sort(a, Opt) {

    // Return empty array
    if (a.length == 0) { return {'data': a, "Ind": new Array(a.length).fill().map((v, i) => i)}; }

    // Decleratyion of variables 
    let i, j, b=[], n, Ind=[], tmp, IA = infoArray(a);

    // Check default input agruments 
    if (Opt == null) { Opt = "ASC"; };

    // Sort array
    if      (isNumber(a))     { return a; }
    else if (isComplexNum(a)) { return a; }
    else if (IA.is1D) {
        b = Copy(a).map((v,i) => {if (isComplexNum(v)) {return Abs(v)} else {return v} } );
        n = b.length;
        Ind = new Array(n).fill().map((v, i) => i);

        // Sort by Ascending order by default
        for (i = 0; i < n; i++) {
            for (j = i + 1; j < n; j++) {
                if (b[j] < b[i]) {
                    tmp = b[i];
                    b[i] = b[j];
                    b[j] = tmp;

                    tmp = Ind[i];
                    Ind[i] = Ind[j];
                    Ind[j] = tmp;
                }
            }
        }
        if      (Opt.toUpperCase() == "ASC")  { return {'data': GetRange(a, Ind),       'Ind': Ind       } }
        else if (Opt.toUpperCase() == "DESC") { return {'data': Flip(GetRange(a, Ind)), 'Ind': Flip(Ind) } }
    }
    else if (IA.is2D) { return a.map((v,i) => Sort(v, Opt) ); }
}
function Sum(a) {
   // Returns the sum of an array
   // If 2D, returns the sum of each row.
   let i, t, IA = infoArray(a);
   if ((isNumber(a)) || (isComplexNum(a))) { return a; }
   else if (IA.is1D) { t=0; for (i=0; i<a.length; i++) { t=Add(a[i],t) }; return t; }
   else if (IA.is2D) { return a.map((v,i) => Sum(v)); }
}
function SwapRowCol(a, rInd1, rInd2, cInd1, cInd2) {
   // Swaps columns and rows at at the same time (in-place swap)
   // rInd1 and rInd2 are arrays containing the row index numbers to be swapped
   // cInd1 and cInd2 are arrays containing the column index numbers to be swapped
   return SwapCol(SwapRow(a, rInd1, rInd2 ), cInd1, cInd2);
}
function SwapRow(a, Ind1, Ind2) {
   // Swaps rows of Array-a[] (in-place swap)
   // Ind1 and Ind2 are arrays containing the row index numbers to be swapped
   let i, temp=[], IA = infoArray(a), res=Copy(a);
   if ((IA.is1D) || (IA.is2D)) {
       for (i=0; i<Ind1.length; i++) {
           temp           = res[ Ind1[i] ];
           res[ Ind1[i] ] = res[ Ind2[i] ];
           res[ Ind2[i] ] = temp;
       }
   }
   return res;
}
function SwapCol(a, Ind1, Ind2) {
   // Swaps columns of Array-a[] (in-place swap)
   // Ind1 and Ind2 are arrays containing the column index numbers to be swapped
   let i, IA = infoArray(a);
   if      (IA.is1D) { return SwapRow(a, Ind1, Ind2); }
   else if (IA.is2D) { for (i=0; i<a.length; i++) { a[i] = SwapRow(a[i], Ind1, Ind2); }}
   return a;
}
function CumSum(a) {
   let i, cs, IA = infoArray(a);
   if ((isNumber(a)) || (isComplexNum(a))) { return a; }
   else if (IA.is1D) { cs=0; return a.map((v,i) => { cs = Add(cs, v); return cs; }) }
   else if (IA.is2D) { return a.map((v,i) => CumSum(v)); }
}
function Sqrt(a) {
   let i, IA = infoArray(a);
   if      (isNumber(a))     { return Math.sqrt(a); }
   else if (isComplexNum(a)) { return a.Sqrt();     }
   else if (IA.is1D)         { return a.map((v,i) => Sqrt(v) ); }
   else if (IA.is2D)         { return a.map((v,i) => Sqrt(v) ); }
}
function Tan(a) {
   let i, IA = infoArray(a);
   if      (isComplexNum(a)) { return a.Tan(); }
   else if (isNumber(a))     { return Math.tan(a); }
   else if (IA.is1D)         { return a.map((v,i) => Tan(v) ); }
   else if (IA.is2D)         { return a.map((v,i) => Tan(v) ); }
}
function Tanh(a) {
   let i, IA = infoArray(a);
   if      (isComplexNum(a)) { return a.Tanh(); }
   else if (isNumber(a))     { return Math.tanh(a); }
   else if (IA.is1D)         { return a.map((v,i) => Tanh(v) ); }
   else if (IA.is2D)         { return a.map((v,i) => Tanh(v) ); }
}
function Transpose(a) {
   // Returns the transpose of A[,] matrix
   // This function is double checked for correctness 2023-06-09
   let IA = infoArray(a);
   if      (isNumber(a))     { return a; }
   else if (isComplexNum(a)) { return a.Conj(); }
   else if (IA.is1D)         { return a.map((v,i) => [Transpose(v)] ); }
   else if (IA.is2D)         { return a[0].map((v,i) =>  a.map((vv,j) =>  Transpose(a[j][i]) )); }
}
function Triu(a) {
   // Returns the upper triangular of Array-a[]
   let i, IA = infoArray(a);
   if (isNumber(a) || isComplexNum(a)) { return a; }
   if (IA.is1D) { return a; }
   else if (IA.is2D) { return a.map((vv,i) => a[i].map((v,j) =>  {if (j>=i) {return v;} else {return 0;}} )); }
}
function Truncate(a, N, Opt) {
   // If Opt is true, first N elements of the data[] array is returned.
   // if Opt is false, last N elements of the data[] array is returned.
   let IA = infoArray(a);
   if (Opt == null) { Opt = true; }
   if      (IA.is1D) { if ((N >= 0) && (N < a.length)) { if (Opt) { return a.slice(0, N); }  else { return a.slice(-N,); } } else { return a; } }
   else if (IA.is2D) { return a.map((v,i) => Truncate(v, N, Opt)); }
}
function Var(a) {
   let IA = infoArray(a);
   if      (isNumber(a) || isComplexNum(a)) { return 0; }
   else if (IA.is1D)                        { return Divide(Sum(Pow(Abs(Subtract(a,Mean(a))),2)),a.length-1); }
   else if (IA.is2D)                        { return a.map((v,i) => Var(v)); }
}
function WilkinsonShift(a, b, c) {
   // Calculate Wilkinson's shift for symmetric matrices:
   // Wilkinson shift is defined as the eigenvalue of the the lower rightmost 2×2 sub-matrix A

   var del, mu, si;

   del = (a - b)/2;
   if (del != 0) { si = Math.sign(del); } else { si = 1;}
   mu = b - si * c * c / (Math.abs(del) + Math.sqrt( del*del + c*c ))

   return mu;
}
function Zeros(a, b) {
   if      (b == null) { return new Array(a).fill(0); }
   else if (b != null) { return new Array(a).fill().map(() => Zeros(b)); }
}
function ZeroPad(a, N) {
   // Zero-Padding to return a new array of length N
   // N must be greater than the length of data[] array; otherwise, the data[] array will be return with no padding.
   let df, IA = infoArray(a);
   if (N == null) { N = a.length; }
   if      (IA.is1D) { df = N - a.length; if (df > 0) { return Copy(a).concat(new Array(df).fill(0)); } else { return a; } }
   else if (IA.is2D) { return a.map((v,i) => ZeroPad(v, N) ); }
}

class ComplexNum {

   constructor(Re, Im) {
       if ((Re != null) && (Im == null)) {
           this.Re = undefined;
           this.Im = undefined;
       }
       else {
           this.Re = Re;  // The real part of the number
           this.Im = Im;  // The imaginary part of the number
       }
   }

   Abs() { return Math.sqrt(Math.pow(this.Re, 2) + Math.pow(this.Im, 2)); };

   Add(a) {
       if      (isNumber(a))     { return new ComplexNum(this.Re + a,    this.Im       ); }
       else if (isComplexNum(a)) { return new ComplexNum(this.Re + a.Re, this.Im + a.Im); }
   };

   Angle() {
       // Return four quadrant arc-tangent of a complex number
       if       (this.Re > 0)                      { return Math.atan(this.Im / this.Re);           }
       else if ((this.Re <  0) && (this.Im >= 0))  { return Math.atan(this.Im / this.Re) + Math.PI; }
       else if ((this.Re <  0) && (this.Im < 0))   { return Math.atan(this.Im / this.Re) - Math.PI; }
       else if ((this.Re == 0) && (this.Im > 0))   { return Math.PI / 2;                            }
       else if ((this.Re == 0) && (this.Im < 0))   { return -Math.PI / 2;                           }
       else if ((this.Re == 0) && (this.Im == 0))  { return 0;                              }
   }

   Conj() { return new ComplexNum(this.Re, -this.Im); }

   Cos() {
       let c1 = Math.cos(this.Re) * Math.cosh(this.Im);
       let c2 = Math.sin(this.Re) * Math.sinh(this.Im);
       return new ComplexNum(c1, -c2);
   }

   Cosh() {
       let c1 = Math.cosh(this.Re) * Math.cos(this.Im);
       let c2 = Math.sinh(this.Re) * Math.sin(this.Im);
       return new ComplexNum(c1, c2);
   }

   Cot() { return this.Tan().Pow(-1); }

   Coth() { return this.Cosh().Divide(this.Sinh()); }

   Deg2Rad() { return new ComplexNum(this.Re*Math.PI/180, this.Im*Math.PI/180 ); }

   Divide(a) {
       if      (isNumber(a))     { return new ComplexNum(this.Re / a, this.Im / a); }
       else if (isComplexNum(a)) {
           let Bj = a.Conj();
           a = a.Multiply(Bj);
           let Nu = this.Multiply(Bj)
           return new ComplexNum(Nu.Re / a.Re, Nu.Im / a.Re);
       }
   }

   Exp() { return new ComplexNum(Math.cos(this.Im), Math.sin(this.Im)).Multiply(Math.exp(this.Re)); }

   IsNaN() { if (Number.isNaN(this.Re) || Number.isNaN(this.Im)) { return true; } else { return false; } }

   IsEqual(a) {
       if      (isNumber(a)) { return this.IsEqual(new ComplexNum(a, 0)); }
       else if (isComplexNum(a)) {
           if ((this.Re==a.Re) && (this.Im==a.Im)) { return true; } else { return false; }
       }
   }

   Imag() { return this.Im; }

   Log() { return new ComplexNum(Math.log(this.Abs()), Math.atan2(this.Im, this.Re)); }

   Log10() { return new ComplexNum(Math.log(this.Abs()) / Math.log(10), Math.atan2(this.Im, this.Re) / Math.log(10)); }

   Multiply(a) {
       if      (isNumber(a))     { return new ComplexNum(this.Re * a, this.Im * a); }
       else if (isComplexNum(a)) { return new ComplexNum(this.Re * a.Re - this.Im * a.Im, this.Re * a.Im + this.Im * a.Re); }
   };

   Pow(a) { return this.Log().Multiply(a).Exp(); }

   Real() { return this.Re; }

   Sec() { return this.Cos().Pow(-1); }

   Sign() { return this.Divide(this.Abs()) }

   Sin() { return new ComplexNum(Math.sin(this.Re) * Math.cosh(this.Im), Math.cos(this.Re) * Math.sinh(this.Im)); }

   Sinh() { return new ComplexNum(Math.sinh(this.Re) * Math.cos(this.Im), Math.cosh(this.Re) * Math.sin(this.Im)); }

   Sqrt() { return this.Pow(0.5); }

   Subtract(a) {
       if      (isNumber(a))     { return new ComplexNum(this.Re - a, this.Im) }
       else if (isComplexNum(a)) { return new ComplexNum(this.Re - a.Re, this.Im - a.Im); }
   };

   Tan() {
       let c1 = new ComplexNum(Math.sin(2 * this.Re), Math.sinh(2 * this.Im));
       let c2 = Math.cos(2 * this.Re) + Math.cosh(2 * this.Im);
       return c1.Divide(c2);
   }

   Tanh() {
        let c1 = new ComplexNum(Math.tanh(this.Re), Math.tan(this.Im));
        let c2 = new ComplexNum(1, Math.tanh(this.Re) * Math.tan(this.Im));
        return c1.Divide(c2);
   }

}
