import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np
import tkinter as tk
from Bio.Seq import Seq
from Bio import pairwise2
from PIL import Image, ImageTk
import io

def fig2img(fig):
    buf = io.BytesIO()
    fig.savefig(buf)
    buf.seek(0)
    img = Image.open(buf)
    return img

def match_score(alpha, beta,match_award,mismatch_penalty):
    if alpha == beta:
        return match_award
    else:
        return mismatch_penalty

def needleman_wunsch(seq1, seq2 ,match_award,mismatch_penalty,gap_penalty):
    m, n = len(seq2), len(seq1)  # length of two sequences
    dt = np.dtype([('diagonal', np.str, 1),
                   ('up', np.str, 1), ('left', np.str, 1)])

    pt_mat = np.zeros((m + 1, n + 1), dtype=dt)

    score_matrix = np.zeros((m + 1, n + 1), dtype=int)     # the DP table
    # Calculate DynamicProgramming table
    for i in range(0, m + 1):
        score_matrix[i][0] = gap_penalty * i
        pt_mat[i][0]['up'] = 'U'
    for j in range(0, n + 1):
        score_matrix[0][j] = gap_penalty * j
        pt_mat[0][j]['left'] = 'L'
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            diagonal = score_matrix[i - 1][j - 1] + \
                match_score(seq2[i - 1], seq1[j - 1],match_award,mismatch_penalty)
            up = score_matrix[i - 1][j] + gap_penalty
            left = score_matrix[i][j - 1] + gap_penalty
            max_pointer = max(diagonal, up, left)
            score_matrix[i][j] = max_pointer
            if (diagonal == max_pointer):
                pt_mat[i][j]['diagonal'] = 'D'
            if (up == max_pointer):
                pt_mat[i][j]['up'] = 'U'
            if (left == max_pointer):
                pt_mat[i][j]['left'] = 'L'

    # Traceback and compute the alignment
    align1, align2 = '', ''
    i, j = m, n  # start from the bottom right cell
    a = np.array([[0, 0, 0, 0]])

    first = True
    global arrows268435459115230927
    while i > 0 and j > 0:  # end toching the top or the left edge
        score_current = score_matrix[i][j]
        score_diagonal = score_matrix[i - 1][j - 1]
        score_up = score_matrix[i][j - 1]
        score_left = score_matrix[i - 1][j]
        #print("current ",i,j,end='')
        a[:, :2] = j, i
        if score_current == score_diagonal + match_score(seq2[i - 1], seq1[j - 1],match_award,mismatch_penalty):
            align1 += seq2[i - 1]
            align2 += seq1[j - 1]
            i -= 1
            j -= 1
        elif score_current == score_left + gap_penalty:
            align1 += seq2[i - 1]
            align2 += '-'
            i -= 1
        elif score_current == score_up + gap_penalty:
            align1 += '-'
            align2 += seq1[j - 1]
            j -= 1
        a[:, 2:] = j, i
        if first:
            # First loop copy a
            arrows = np.copy(a)
            first = False
        else:
            # Concatenate origin-target
            arrows = np.concatenate((arrows, a), axis=0)

    a[:, :2] = a[:, 2:]

    # Finish tracing up to the top left cell
    while i > 0:
        align1 += seq2[i - 1]
        align2 += '-'
        i -= 1
        a[:, 2:] = j, i
        arrows = np.concatenate((arrows, a), axis=0)
    while j > 0:
        align1 += '-'
        align2 += seq1[j - 1]
        j -= 1
        a[:, 2:] = j, i
        arrows = np.concatenate((arrows, a), axis=0)
    return score_matrix, pt_mat

def getoutput(seq1,seq2,match_award,mismatch_penalty,gap_penalty):
    seq1,seq2 = sorted([seq1,seq2],key=len)[::-1]
    plt.rcParams["figure.figsize"] = 4, 5
    param = {"grid.linewidth": 1.6,
             "grid.color": "lightgray",
             "axes.linewidth": 1.6,
             "axes.edgecolor": "lightgray"}
    plt.rcParams.update(param)
    print('here')

    # Data
    headh = seq1
    headv = seq2

    v, pt_mat = needleman_wunsch(seq1,seq2,match_award,mismatch_penalty,gap_penalty)
    print('here2')

    # Plot
    fig, ax = plt.subplots()
    ax.set_xlim(-1.5, v.shape[1] - .5)
    ax.set_ylim(-1.5, v.shape[0] - .5)
    ax.invert_yaxis()
    for i in range(v.shape[0]):
        for j in range(v.shape[1]):
            ax.text(j, i, v[i, j], ha="center", va="center")
    for i, l in enumerate(headh):
        ax.text(i + 1, -1, l, ha="center", va="center", fontweight="semibold")
    for i, l in enumerate(headv):
        ax.text(-1, i + 1, l, ha="center", va="center", fontweight="semibold")

    ax.xaxis.set_minor_locator(ticker.FixedLocator(
        np.arange(-1.5, v.shape[1] - .5, 1)))
    ax.yaxis.set_minor_locator(ticker.FixedLocator(
        np.arange(-1.5, v.shape[1] - .5, 1)))
    plt.tick_params(axis='both', which='both', bottom='off', top='off',
                    left="off", right="off", labelbottom='off', labelleft='off')
    ax.grid(True, which='minor')


    arrowprops = dict(facecolor='black', alpha=0.5, lw=0,
                      shrink=0.2, width=2, headwidth=5, headlength=7)
    print('here3')
    # all paths
    for i in range(1, pt_mat.shape[0]):
        for j in range(1, pt_mat.shape[1]):
            if(pt_mat[i][j]['left'] != ''):
                ax.annotate("", xy=(j - 1, i),
                            xytext=(j, i), arrowprops=arrowprops)
            if(pt_mat[i][j]['diagonal'] != ''):
                ax.annotate("", xy=(j - 1, i - 1),
                            xytext=(j, i), arrowprops=arrowprops)
            if(pt_mat[i][j]['up'] != ''):
                ax.annotate("", xy=(j, i - 1),
                            xytext=(j, i), arrowprops=arrowprops)
    print('here4')
    # optimal path
    arrowprops.update(facecolor='cyan')
    for i in range(arrows.shape[0]):
        ax.annotate("", xy=arrows[i, 2:],
                    xytext=arrows[i, :2], arrowprops=arrowprops)

    im = fig2img(fig)
    alignments = ''
    for a in pairwise2.align.globalms(seq1, seq2, match_award, mismatch_penalty, gap_penalty, 0):
                alignments+=pairwise2.format_alignment(*a, full_sequences=True)+'\n'
    return im, alignments
    

arrows = None
im = None


master = tk.Tk()

# TK VARS
boxvar = tk.StringVar()
boxvar2 = tk.StringVar()
boxvar3 = tk.StringVar()
boxvar4 = tk.StringVar()
boxvar5 = tk.StringVar()

# BUTTON FUNCTIONS
def showresults():
    im, alignments = getoutput(box.get(),box2.get(),int(box3.get()),int(box4.get()),int(box5.get()))
    img = ImageTk.PhotoImage(im.resize([400,400],resample=1))
    canvas.create_image(0,0,anchor=tk.NW, image=img)
    canvas.image = img
    txt.delete('1.0', tk.END)
    txt.insert('1.0', alignments)
    boxvar.set('')
    boxvar2.set('')


# UI - INPUT MENU
inputmenu = tk.Frame(master)
label = tk.Label(inputmenu, text="Sequence 1")
label2 = tk.Label(inputmenu, text="Sequence 2")
label3= tk.Label(inputmenu, text="Match")
label4 = tk.Label(inputmenu, text="Mismatch")
label5 = tk.Label(inputmenu, text="Gap")
box = tk.Entry(inputmenu, textvariable = boxvar)
box2 = tk.Entry(inputmenu, textvariable = boxvar2)
box3 = tk.Entry(inputmenu, textvariable = boxvar3)
box4 = tk.Entry(inputmenu, textvariable = boxvar4)
box5 = tk.Entry(inputmenu, textvariable = boxvar5)
btn = tk.Button(inputmenu, text = 'Display\nResults', command=showresults)

# UI - TEXT OUTPUT
textoutput = tk.Frame(master)
txt = tk.Text(textoutput, width=30)
txt.delete('1.0', tk.END)
txt.insert('1.0', '\n'*10+'\tNo Input!')
scrollb = tk.Scrollbar(textoutput, command=txt.yview)
txt['yscrollcommand'] = scrollb.set

# UI - IMAGE OUTPUT
canvas = tk.Canvas(master, height=400)
imd = Image.fromarray(np.random.randint(low=100,high=200,size=(400,400)),'I')
img = ImageTk.PhotoImage(imd)
canvas.create_image(0,0,anchor=tk.NW, image=img)
canvas.image = img

# GRIDDING - INPUT MENU
label.grid(row=0,column=0)
label2.grid(row=0,column=1)
label3.grid(row=0,column=2)
label4.grid(row=0,column=3)
label5.grid(row=0,column=4)
box.grid(row=1,column=0)
box2.grid(row=1,column=1)
box3.grid(row=1,column=2)
box4.grid(row=1,column=3)
box5.grid(row=1,column=4)
btn.grid(row=0, rowspan=2, column=5)
inputmenu.pack(side='top')

# GRIDDING - TEXT OUTPUT
txt.grid(row=0, column=0, sticky='nsew', padx=2, pady=2)
scrollb.grid(row=0, column=1, sticky='nsew')
textoutput.pack(side='left')

# GRIDDING - IMAGE OUTPUT
canvas.pack(side='right')

master.mainloop()

