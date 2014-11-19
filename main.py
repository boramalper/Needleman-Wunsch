# Copyright (c) 2014 Mert Bora Alper <boramalper@gmail.com>
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.

import copy
from configuration import *

def init():
    similarity_matrix[1][0] = similarity_matrix[0][1]  # AG = GA
    similarity_matrix[2][0] = similarity_matrix[0][2]  # AC = CA
    similarity_matrix[2][1] = similarity_matrix[1][2]  # GC = CG
    similarity_matrix[3][0] = similarity_matrix[0][3]  # AC = CA
    similarity_matrix[3][1] = similarity_matrix[1][3]  # AC = CA
    similarity_matrix[3][2] = similarity_matrix[2][3]  # AC = CA

    global sm
    sm = {}

    for yi, y in enumerate("AGCT"):
        sm[y] = {}
        for xi, x in enumerate("AGCT"):
            sm[y][x] = similarity_matrix[yi][xi]


def main():
    init()

    dna1 = input("DNA1< ").upper()
    dna2 = input("DNA2< ").upper()

    table = create_table(dna1, dna2)
    print_table(dna1, dna2, table)

    matches = find_alignments(dna1, dna2, table)
    print_alignments(matches)


def create_table(dna1, dna2):
    # create table
    table = []
    for y in range(len(dna2) + 1):
        table.append([])
        for x in range(len(dna1) + 1):
            table[y].append({"value": 0, "arrows": [False, False, False]})
            #               ^ A cell of the Table.   Left   Diag.  Up

    # fill first row
    for i in range(len(dna1) + 1):
        table[0][i]["value"] = -i
        table[0][i]["arrows"][0] = True

    # fill first column
    for i in range(len(dna2) + 1):
        table[i][0]["value"] = -i
        table[i][0]["arrows"][2] = True

    # while filling first row and column, we allowed this mistake for the sake
    # of simplicity but we need to fix that
    table[0][0]["arrows"] = [False, False, False]

    # calculate the rest
    for i2, b2 in enumerate(dna2, start=1):
        for i1, b1 in enumerate(dna1, start=1):
            left_score = table[i2][i1-1]["value"] + indel
            dia_score = table[i2-1][i1-1]["value"] + sm[b2][b1]
            up_score = table[i2-1][i1]["value"] + indel

            max_score = max(left_score, dia_score, up_score)
            table[i2][i1]["value"] = max_score

            if left_score == max_score:
                table[i2][i1]["arrows"][0] = True
            if dia_score == max_score:
                table[i2][i1]["arrows"][1] = True
            if up_score == max_score:
                table[i2][i1]["arrows"][2] = True

    return table


def print_table(dna1, dna2, table):
        # Numbers Table
    print("Scores:")
    print("\t      ", end="")
    for i in dna1:
        print(i, end="  ")
    print()

    for iy, y in enumerate(table):
        print("\t%s" % (" " + dna2)[iy], end=" ")
        for x in y:
            print("% d" % x["value"], end=" ")
        print()

    print()

    # Arrows Table
    print("Arrows: (value = 1 if left + 2 if diagonal + 4 if up)")
    print("\t      ", end="")
    for i in dna1:
        print(i, end="  ")
    print()

    for iy, y in enumerate(table):
        print("\t%s" % (" " + dna2)[iy], end=" ")
        for x in y:
            summe = x["arrows"][0] + x["arrows"][1]*2 + x["arrows"][2]*4
            print("% d" % summe, end=" ")
        print()

    print("\n")


def find_alignments(dna1, dna2, table):
    alignments = []

    # trace_r: trace arrows back to their origin
    # !recursive function!
    def trace_r(y, x, alignment):
        if table[y-1][x-1]["arrows"][0]:  # LEFT
            alignment_c = copy.deepcopy(alignment)  # create a branch
            alignment_c["DNAs"][0] += dna1[x-1]
            alignment_c["DNAs"][1] += "-"
            alignment_c["score"] += table[y-1][x-1]["value"]
            trace_r(y, x-1, alignment_c)
            del alignment_c

        if table[y-1][x-1]["arrows"][1]:  # DIAGONAL
            alignment_c = copy.deepcopy(alignment)  # create a branch
            alignment_c["DNAs"][0] += dna1[x-1]
            alignment_c["DNAs"][1] += dna2[y-1]
            alignment_c["score"] += table[y-1][x-1]["value"]
            trace_r(y-1, x-1, alignment_c)
            del alignment_c

        if table[y-1][x-1]["arrows"][2]:  # UP
            alignment_c = copy.deepcopy(alignment)  # create a branch
            alignment_c["DNAs"][0] += "-"
            alignment_c["DNAs"][1] += dna2[y-1]
            alignment_c["score"] += table[y-1][x-1]["value"]
            trace_r(y-1, x, alignment_c)
            del alignment_c

        # origin :)
        if y-1 == -len(table) and x-1 == -len(table[0]):
            alignment["DNAs"] = alignment["DNAs"][0][::-1], alignment["DNAs"][1][::-1]
            score = calculate_score(alignment["DNAs"][0], alignment["DNAs"][1], alignment["score"])
            alignments.append({"DNAs": alignment["DNAs"], "score": score})

    trace_r(0, 0, {"DNAs": ["", ""], "score": 0})

    alignments = sorted(alignments, key=lambda match: -match["score"])

    return alignments


def print_alignments(alignments):
    print("%d Alignments:" % (len(alignments)))
    for a in alignments:
        print("\t%s\tScore:" % (a["DNAs"][0]))
        print("\t%s\t%d" % (a["DNAs"][1], a["score"]))
        print()


def calculate_score(alg1, alg2, score):
    for i in range(len(alg1)):
        if alg1[i] == '-':
            if i >= 1 and alg1[i-1] == '-':
                score += (extend_indel - indel)
            else:
                score += (new_indel - indel)

    for i in range(len(alg2)):
        if alg2[i] == '-':
            if i >= 1 and alg2[i-1] == '-':
                score += (extend_indel - indel)
            else:
                score += (new_indel - indel)

    return score

if __name__ == "__main__":
    main()
