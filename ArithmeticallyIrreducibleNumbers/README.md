A digit string S in base b is a finite string made of digits from 0, 1, ..., b-1. We are interested in some arithmetic properties P a digit string may satisfy. Here are some, with example / non-example digit strings in base 10.

| Name       | Description                               | Example           | NonExample |
| ---------- | ----------------------------------------- | ----------------- | ---------- |
| SumBalance | Can be split into two equal sums.         | 74524: 7+4=5+2+4  | 39461      |
| AS         | Some expression using +- equals zero.     | 9461: 9-4-6+1=0   | 1248       |
| ASMDP      | Some expression using +-x/() equals zero. | 2958: 2x(9-5)-8=0 | 8985898    |

(Note: we always define the string 0 to have the property.) For each of the above P, there are inifinitely many examples. For example, 11, 1111, 111111, or any even amount of 1's. So the strings satisying P are not too interesting to me, but something else is.

Given such a property P, we say that a digit string S is P-irreducible if no contiguous substring of S satisfies P. For each above P, except maybe the first, there are only finitely many P-irreducible digit strings! We strive to count them.

This is well-suited to a back-tracking approach. If there are finitely many P irreducible strings, there is a maximum length, which we denote Lmax. At each depth of the backtracking algorithm, there are at most b children states, and the depth is at most Lmax, so the stack has at most b\*Lmax entries at any time.

Under this approach, we have verified some string abcde is irreducible, and we ask for what x is abcdex irreducible. Notice we do not need to check every substring of abcdex, as most have already been checked. We need only check the suffixes x, ex, dex, cdex, bcdex, abcdex. If P is nice, we can usually check all suffixes in the same calculation, maybe even all candidates digits x too.

Here are some known results:

<details>
<summary>SumBalance-irreducible</summary>

| base | total   | length | time(s) |
| ---- | ------- | ------ | ------- |
| 2    | 1       | 1      | 0.000   |
| 3    | 6       | 3      | 0.000   |
| 4    | 33      | 7      | 0.001   |
| 5    | 310     | 9      | 0.001   |
| 6    | 5188    | 15     | 0.052   |
| 7    | 257385  | 21     | 1.978   |
| 8    | 9699463 | 31     | 82.699  |

Questions:

- Can we prove this is always finite / bound the max length?
- Can we calculate terms 9 and 10?

Term 9 takes over 2 hours.

</details>

<details>
<summary>AS-irreducible</summary>

| base | total    | length | time(s) |
| ---- | -------- | ------ | ------- |
| 2    | 1        | 1      | 0.000   |
| 3    | 5        | 3      | 0.000   |
| 4    | 14       | 3      | 0.000   |
| 5    | 84       | 7      | 0.001   |
| 6    | 210      | 7      | 0.001   |
| 7    | 993      | 7      | 0.003   |
| 8    | 2310     | 7      | 0.006   |
| 9    | 32318    | 15     | 0.279   |
| 10   | 48469    | 15     | 0.304   |
| 11   | 269405   | 15     | 1.887   |
| 12   | 396146   | 15     | 2.052   |
| 13   | 6368008  | 15     | 38.143  |
| 14   | 7972901  | 15     | 42.715  |
| 15   | 41304188 | 15     | 230.374 |
| 16   | 51298219 | 15     | 274.692 |

Proposition: The maximum length is $<2b$.

Questions:

- Can we calculate more terms?
- Why does the maximum length seem to be given by [A003817](https://oeis.org/A003817)?
- Why the jump in total and runtime on odd bases?
</details>

I have some other thoughts but this will do for now.
