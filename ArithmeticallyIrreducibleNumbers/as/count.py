#! /usr/bin/env python
"""
Back-tracking algorithms to count the AS-irreducible strings in base b.
We know L_max < 2*b.
"""

def count_AS_irreducible_1(base):
    """
    O(b L_max^2) time per node, O(b L_max) space.
    Calculates candidates from scratch each node.
    """
    digits = tuple(range(1, base))
    counts = [0] * (2 * base)
    counts[0] = base-1  # all the single digits
    arr = []
    tm1 = -1    # total minus 1
    blank_mask = (1<<base)-2    # all non-zero digits
    # stack will have at most b*L_max states

    stack = [(x, False) for x in digits]
    while stack:
        digit, second_pass = stack.pop()

        # backtrack changes
        if second_pass:
            arr.pop()
            tm1 -= digit
            continue
        stack.append((digit, True))

        # enact changes
        arr.append(digit)
        tm1 += digit

        # find the digits x such that no suffix of arr+[x] has the AS property,
        # ie. exists signed summation sigma such that sigma(suffix)=0.
        # O(bL) space, O(bL^2) time using bit manipulations
        rz = tm1 // 2   # right zeros
        sums = 1 << rz
        cm = blank_mask << rz
        for y in reversed(arr):
            sums = (sums << y) | (sums >> y)
            cm &= ~sums
            if not cm:
                break
        else:
            # extract the valid digits
            cm >>= rz
            counts[len(arr)] += cm.bit_count()
            while cm:
                lsb = cm & -cm
                x = lsb.bit_length() - 1
                stack.append((x, False))
                cm ^= lsb

    # remove the unused lengths
    while not counts[-1]:
        counts.pop()
    return counts


def count_AS_irreducible_2(base):
    """
    O(b L_max) time per node, O(b L_max^2) space.
    Stores a list of bitmasks to track the disallowed candidates.
    """
    digits = tuple(range(1, base))
    counts = [0] * (2 * base)
    counts[0] = base-1  # all the single digits
    arr = []
    rz = (base-1) * base # right zeros, <=(b-1) L_max // 2 < (b-1)b
    zero = 1 << rz
    sums_list = [zero]
    blank_mask = (1<<base)-2    # all non-zero digits
 
    # stack will have at most b*L_max states
    stack = [(x, False) for x in digits]
    while stack:
        digit, second_pass = stack.pop()

        # backtrack changes
        if second_pass:
            arr.pop()
            sums_list.pop()
            continue
        stack.append((digit, True))

        # enact changes
        arr.append(digit)
        s = sums_list[-1]
        s = (s << digit) | zero | (s >> digit)
        sums_list.append(s)

        # extract the valid digits
        cm = blank_mask & ~ (s >> rz)
        counts[len(arr)] += cm.bit_count()
        while cm:
            lsb = cm & -cm
            x = lsb.bit_length() - 1
            stack.append((x, False))
            cm ^= lsb

    # remove the unused lengths
    while not counts[-1]:
        counts.pop()
    return counts


if __name__ == '__main__':
    import time
    print(f"base        total  length     time(s)")
    print("-------------------------------------")
    for base in range(2, 17):
        st = time.time()
        counts = count_AS_irreducible_2(base)
        print(f"{base:>4} {sum(counts):>12}   {len(counts):>5}  {time.time()-st:>10.3f}")
