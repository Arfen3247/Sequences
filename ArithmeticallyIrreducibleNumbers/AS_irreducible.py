"""
A back-tracking algorithm to count the AS-irreducible strings in base b.
Uses O(b L_max) space. We know L_max < 2*b.
"""


BIT_POS = [tuple(i for i in range(8) if (byte >> i) & 1) for byte in range(256)]
def count_AS_irreducible(base):
    digits = tuple(range(1, base))
    counts = [0] * (2 * base)
    counts[1] = base-1  # all the single digits
    arr = []
    length = 1

    bounds = [base * (length+1)//2 for length in range(2*base)]
    blank_mask = (1<<base)-2

    # stack will have at most b*L_max states
    stack = [(x, False) for x in digits]
    while stack:
        digit, second_pass = stack.pop()

        # backtrack changes
        if second_pass:
            length -= 1
            arr.pop()
            continue
        stack.append((digit, True))

        # enact changes
        length += 1
        arr.append(digit)

        # find the digits x such that no suffix of arr+[x] has the AS property,
        # ie. exists signed summation sigma such that sigma(suffix)=0.
        # O(bL) space, O(bL^2) time using bit manipulations (need 2*b*L bits, maybe b*L).
        cm = blank_mask
        bound = bounds[length]
        sums = 1 << bound
        for y in reversed(arr):
            sums = (sums << y) | (sums >> y)
            cm &= ~(sums>>bound)
            if not cm:
                break
        if not cm:
            continue

        # read off the allowed candidates x in O(b) time, add to stack
        # this part is probably dumb
        offset = 0
        for byte_val in cm.to_bytes((cm.bit_length() + 7) // 8, 'little'):
            for bit in BIT_POS[byte_val]:
                stack.append((offset+bit, False))
            offset += 8
        counts[length] += cm.bit_count()

    # remove the unused lengths
    length = len(counts) - 1
    while not counts[length]:
        length -= 1
    return counts

import time
if __name__ == '__main__':
    print(f"base        total  length     time(s)")
    print("-------------------------------------")
    for base in range(2, 17):
        st = time.time()
        counts = count_AS_irreducible(base)
        print(f"{base:>4} {sum(counts):>12}   {len(counts)-1:>5}  {time.time()-st:>10.3f}")

"""
base        total  length     time(s)
-------------------------------------
   2            1       1       0.000
   3            5       3       0.000
   4           14       3       0.000
   5           84       7       0.001
   6          210       7       0.001
   7          993       7       0.003
   8         2310       7       0.006
   9        32318      15       0.279
  10        48469      15       0.304
  11       269405      15       1.887
  12       396146      15       2.052
  13      6368008      15      38.143
  14      7972901      15      42.715
  15     41304188      15     230.374
  16     51298219      15     274.692
"""