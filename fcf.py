import math as m
import matplotlib.pyplot as plt
import numpy as np
import scipy.constants as c
from scipy.special import assoc_laguerre
import tkinter as tk
from tkinter import *
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
import matplotlib.figure as f
import periodictable

### Reference DOI: 10.1063/1.443949

# Quantum Numbers: v = v', V = v''

def calculateWBar(w_gs, w_es):
    w_bar = 4 * (w_gs * w_es) / (w_gs + 2 * m.sqrt(w_gs * w_es) + w_es)
    return w_bar

def q(v, V, u):
  z = (m.factorial(v) / m.factorial(V))
  if (V - v) > -1:
    return z * (u**(V - v)) * m.exp(-u) * (assoc_laguerre(u, v, V - v)**2)
  else:
    return "v_es - v_gs > -1 must be satisfied"

def progression(v, n_V, u):
   arr = []
   for i in range(n_V+1):
        if (i < v):
            arr.append(int(100*q(i, v, u)))
        else:
           arr.append(int(100*q(v, i, u)))
   return arr 

# Graphing

def generateColors(v_gs,color):
    arr = []
    j = 0
    for i in range(v_gs+1):
        if j > 3:
            j = 0
        arr.append(color[j])
        j+=1
    return(arr)

   
def createTable():
    
    secondary_window = tk.Toplevel(bg='white')
    secondary_window.title("Individual Progressions")
    secondary_window.geometry("600x500")

    main_frame = Frame(secondary_window)
    main_frame.pack(fill=BOTH, expand=1)
    
    # Update Constants
    r_gs = e4.get()
    r_es = e3.get()

    delta_r = float(r_es) - float(r_gs)

    w_gs = e6.get()
    w_es = e5.get()

    w_bar = calculateWBar(float(w_gs), float(w_es))

    eA = periodictable.elements[int(e1.get())]
    eB = periodictable.elements[int(e2.get())]

    m1 = eA.mass
    m2 = eB.mass
    
    el1 = eA.symbol
    el2 = eB.symbol

    l1.config(text=el1)
    l2.config(text=el2)

    mol = eA.symbol + "-" + eB.symbol

    mu = (m1 * m2) / (m1 + m2)
    u = (mu * w_bar) * (delta_r)**2 / 67.4425

    v_gs = int(e8.get())
    v_es = int(e7.get())
    
    xs = [i for i in range(v_gs+1)]
    lst = []
    for i in range(v_es+1):
        temp = progression(i, v_gs+1,u)
        temp.insert(0,i)
        lst.append(temp)
        print(lst)
        for j in range(v_gs+1):
            e = Entry(main_frame, width=10, fg='blue',font=('Arial',10,'bold'))
            e.grid(row=i, column=j)
            e.insert(END, lst[i][j])
    print(lst)

def plotBars():

    # Creating new window/frame
    
    secondary_window = tk.Toplevel(bg='white')
    secondary_window.title("Individual Progressions")
    secondary_window.geometry("600x500")

    main_frame = Frame(secondary_window)
    main_frame.pack(fill=BOTH, expand=1)

    ## Data processing

    # Update Constants
    r_gs = e4.get()
    r_es = e3.get()

    delta_r = float(r_es) - float(r_gs)

    w_gs = e6.get()
    w_es = e5.get()

    w_bar = calculateWBar(float(w_gs), float(w_es))
    
    eA = periodictable.elements[int(e1.get())]
    eB = periodictable.elements[int(e2.get())]

    m1 = eA.mass
    m2 = eB.mass
    
    el1 = eA.symbol
    el2 = eB.symbol

    l1.config(text=el1)
    l2.config(text=el2)

    mu = (m1 * m2) / (m1 + m2)
    u = (mu * w_bar) * (delta_r)**2 / 67.4425

    v_gs = int(e8.get())
    v_es = int(e7.get())

    xs = [i for i in range(v_gs+1)]

    if v_es == 0:
       fig,ax = plt.subplots()
       
       if el1 == el2:  
            fig.suptitle('Vibronic Progressions for ${}_2$'.format(el1))
       else:
            fig.suptitle('Vibronic Progressions for {}-{}'.format(el1,el2))
        
       ys = progression(0, v_gs, u)
       
    #    f_poly = np.polyfit(xs, ys, deg=v_gs+1)
    #    x_fit = np.linspace(min(xs), max(xs), 100)
    #    y_fit = np.polyval(f_poly, x_fit)
       
       ax.bar(xs,ys, width=0.5)
       ax.set_title('Excited State: v\'= ' + str(0))
       ax.set_xlabel('Ground State (v\'\')')
       ax.set_ylabel('Franck-Condon Factor (%)')
       
    #    plt.plot(x_fit, y_fit, '--')
    else:
        fig = f.Figure(figsize=(17, 12), dpi=80)

        rows = 3
        cols = 0
        n = (v_es + 1)//rows

        if (v_es + 1) % rows == 0:
            cols = n
        else:
            cols = n+1

        if (v_es + 1) <= rows:
            ax = fig.subplots(v_es+1)
            
            if el1 == el2:  
                fig.suptitle('Vibronic Progression for ${}_2$'.format(el1))
            else:
                fig.suptitle('Vibronic Progression for {}-{}'.format(el1,el2))
                
            fig.tight_layout(pad=3.0)
            for i in range(v_es+1):
                ys = progression(i,v_gs,u)
                
                ax[i].bar(xs, ys, width=0.5)
                ax[i].set_title('Excited State: v\'= ' + str(i))
                ax[i].set_xlabel('Ground State (v\'\')')
                ax[i].set_ylabel('Franck-Condon Factor (%)')
                ax[i].set_facecolor('lightgrey')
                # ax[i].plot(x_fit, y_fit, '--')
                
        else:
            ax = fig.subplots(nrows=rows,ncols=cols)
            
            if el1 == el2:  
                fig.suptitle('Vibronic Progressions for ${}_2$'.format(el1))
            else:
                fig.suptitle('Vibronic Progressions for {}-{}'.format(el1,el2))
                
            fig.tight_layout(pad=3.0)
            m, n = 0, 0

            for i in range(v_es+1):
                if m > rows-1:
                    m = 0
                    n += 1 
                ys = progression(i,v_gs,u)
                ax[m,n].bar(xs, ys, width=0.5)
                ax[m,n].set_title('Excited State: v\'= ' + str(i))
                ax[m,n].set_xlabel('Ground State (v\'\')')
                ax[m,n].set_ylabel('Franck-Condon Factor (%)')
                ax[m,n].set_facecolor('lightgrey')
                m += 1

    # Placing 2D Bar Chart

    canvas = FigureCanvasTkAgg(fig, master = secondary_window)

    canvas.get_tk_widget().place(x=0, y=50)

    toolbar = NavigationToolbar2Tk(canvas, secondary_window, pack_toolbar = False)
    toolbar.update()
    toolbar.place(x=0,y=0)

    frame.pack()
    
    canvas.draw()

def plot():
    # Clear Canvas
    ax.clear()
    # Update Constants
    r_gs = e4.get()
    r_es = e3.get()

    delta_r = float(r_es) - float(r_gs)

    w_gs = e6.get()
    w_es = e5.get()

    w_bar = calculateWBar(float(w_gs), float(w_es))

    eA = periodictable.elements[int(e1.get())]
    eB = periodictable.elements[int(e2.get())]

    m1 = eA.mass
    m2 = eB.mass
    
    el1 = eA.symbol
    el2 = eB.symbol

    l1.config(text=el1)
    l2.config(text=el2)

    mu = (m1 * m2) / (m1 + m2)
    u = (mu * w_bar) * (delta_r)**2 / 67.4425

    v_gs = int(e8.get())
    v_es = int(e7.get())

    # Load plot
    
    color_arr = ['r', 'g', 'b', 'y']

    colors = generateColors(v_es, color_arr)
    yticks = [i for i in range(v_es+1)]
    yticks.reverse()

    for c, k in zip(colors, yticks):
        # Create data set in plane k
        xs = np.arange(v_gs+1)
        ys = progression(k, v_gs, u)

        # Coloring each set
        cs = [c] * len(xs)

        # Plot the bar graph given by xs and ys on the plane y=k with 80% opacity.
        ax.bar(xs, ys, zs=k, zdir='y', color=cs, alpha=0.8)
        
    if el1 == el2:  
        ax.set_title('Series of Vibronic Progression(s) for ${}_2$'.format(el1))
    else:
        ax.set_title('Series of Vibronic Progression(s) for {}-{}'.format(el1,el2))

    ax.set_xlabel('Ground State (v\'\')')
    ax.set_ylabel('Excited State (v\')')
    ax.set_zlabel('Franck-Condon Factor (%)')
    ax.set_yticks(yticks)
    
    
        
    canvas.draw()

fig = f.Figure(figsize=(4.5, 4), dpi=100, linewidth=5,edgecolor='darkseagreen')
ax = fig.add_subplot(projection='3d')

# Initialization

root = tk.Tk()
root.geometry('1000x700')
root.title("FCF-Calc")
root.configure(background='dimgray')

# App
frame = tk.Frame(root)
label = tk.Label(root, text = "FCF Calculator", font=("Courier",40,"bold","italic"))
label.configure(background='dimgray',fg='chartreuse')
label.pack()

# Entries
offset1 = 0
offset2 = 100 

l = Label(root, text="Parameters", font=("Courier",25,"bold","italic"))
l.place(x=85 + offset1, y=50 + offset2)
l.configure(background='dimgray',fg='chartreuse')

# Masses

l1 = Label(root, text="Z (1)", font=("Courier",12,"bold"))
l1.place(x=95 + offset1, y=120 + offset2)
l1.configure(background='dimgray',fg='chartreuse')

l2 = Label(root, text="Z (2)", font=("Courier",12,"bold"))
l2.place(x=195+ offset1, y=120+ offset2)
l2.configure(background='dimgray',fg='chartreuse')

e1 = Entry(root, width=10, font=("Courier",9,"bold"))
e1.insert(1,23)
e1.place(x=100+ offset1, y=150+ offset2)

e2 = Entry(root, width=10, font=("Courier",9,"bold"))
e2.insert(1,9)
e2.place(x=200+ offset1, y=150+ offset2)

# Internuclear Distances

l3 = Label(root, text="R'e", font=("Courier",12,"bold"))
l3.place(x=95+ offset1, y=170+ offset2)
l3.configure(background='dimgray',fg='chartreuse')

l4 = Label(root, text="R''e", font=("Courier",12,"bold"))
l4.place(x=195+ offset1, y=170+ offset2)
l4.configure(background='dimgray',fg='chartreuse')

e3 = Entry(root, width=10, font=("Courier",9,"bold"))
e3.insert(1,1.703)
e3.place(x=100+ offset1, y=200+ offset2)

e4 = Entry(root, width=10, font=("Courier",9,"bold"))
e4.insert(1,1.775)
e4.place(x=200+ offset1, y=200+ offset2)

# Vibrational Frequencies

l5 = Label(root, text="w'e", font=("Courier",12,"bold"))
l5.place(x=95+ offset1, y=220+ offset2)
l5.configure(background='dimgray',fg='chartreuse')

l6 = Label(root, text="w''e", font=("Courier",12,"bold"))
l6.place(x=195+ offset1, y=220+ offset2)
l6.configure(background='dimgray',fg='chartreuse')

e5 = Entry(root, width=10, font=("Courier",9,"bold"))
e5.insert(1,1620)
e5.place(x=100+ offset1, y=250+ offset2)

e6 = Entry(root, width=10, font=("Courier",9,"bold"))
e6.insert(1,1386)
e6.place(x=200+ offset1, y=250+ offset2)

# Vibrational QNs

l7 = Label(root, text="upper v'", font=("Courier",12,"bold"))
l7.place(x=95+ offset1, y=270+ offset2)
l7.configure(background='dimgray',fg='chartreuse')

l8 = Label(root, text="upper v''", font=("Courier",12,"bold"))
l8.place(x=195+ offset1, y=270+ offset2)
l8.configure(background='dimgray',fg='chartreuse')

e7 = Entry(root, width=10, font=("Courier",9,"bold"))
e7.insert(1,0)
e7.place(x=100+ offset1, y=300+ offset2)

e8 = Entry(root, width=10, font=("Courier",9,"bold"))
e8.insert(1,8)
e8.place(x=200+ offset1, y=300+ offset2)

# Delta r_e controls (second 2D Bar Chart)

# l7 = Label(root, text="delta r_e (1)", font=("Courier",12,"bold"))
# l7.place(x=95+ offset1, y=270+ offset2)
# l7.configure(background='dimgray',fg='chartreuse')

# l8 = Label(root, text="delta r_e (2)''", font=("Courier",12,"bold"))
# l8.place(x=195+ offset1, y=270+ offset2)
# l8.configure(background='dimgray',fg='chartreuse')

# e7 = Entry(root, width=10, font=("Courier",9,"bold"))
# e7.insert(1,0)
# e7.place(x=100+ offset1, y=300+ offset2)

# e8 = Entry(root, width=10, font=("Courier",9,"bold"))
# e8.insert(1,8)
# e8.place(x=200+ offset1, y=300+ offset2)

# Buttons

b1 = Button(root, text= "Plot Vibronic Progression", command = plot, fg = 'chartreuse', bg = 'darkslategray',font = ("Courier",12,"bold"))
b1.place(x=50,y=500)
b2 = Button(root, text= "Plot Individual Progressions", command = plotBars, fg = 'chartreuse', bg = 'darkslategray',font = ("Courier",12,"bold"))
b2.place(x=50,y=460)
b2 = Button(root, text= "View Tables", command = createTable, fg = 'chartreuse', bg = 'darkslategray',font = ("Courier",12,"bold"))
b2.place(x=50,y=540)

canvas = FigureCanvasTkAgg(fig, master = root)

canvas.get_tk_widget().place(x=400, y=150)

toolbar = NavigationToolbar2Tk(canvas, root, pack_toolbar = False)
toolbar.update()
toolbar.place(x=400, y=550)

frame.pack()

root.mainloop()