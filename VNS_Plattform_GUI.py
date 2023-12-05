
global GUI_version; GUI_version = '1.6'

import os
import copy
import numpy as np
import imageio # version 2.10.5 required!
import tkinter as tk
import tkinter.filedialog
import matplotlib as mpl
import matplotlib.pyplot as plt
from PIL import ImageTk, Image
from scipy.linalg import expm, norm
from mpl_toolkits.axes_grid1 import Divider, Size

def process_entry(entry, result):
    try:
        result.configure(state='normal')
        result.delete(0, tk.END)
        result.configure(state='disabled')
        float(entry.get()) # throws error if not a number
        result.configure(state='normal')
        result.insert(0, "%.2f" % float(entry.get()))
        result.configure(state='disabled')
        if result.get().find("<") == -1:
            set_green(result)
        else:
            set_red(result)
    except:
        pass

def set_result(result, text):
    result.configure(state='normal')
    result.delete(0, tk.END)
    result.insert(0, text)
    result.configure(state='disabled')
    if result.get().find("<") == -1:
        set_green(result)
    else:
        set_red(result)

def set_entry(entry, text):
    entry.delete(0, tk.END)
    entry.insert(0, text)
        
def set_green(entry):
    entry.configure(state='normal')
    entry.configure(disabledbackground='#DCFFDC')
    entry.configure(state='disabled')

def set_red(entry):
    entry.configure(state='normal')
    entry.configure(disabledbackground='#FFDCDC')
    entry.configure(state='disabled')

def rot_around(vec, theta, axis):
    return np.dot(expm(np.cross(np.eye(3), axis/norm(axis)*theta)), vec)

def points_to_coordinates(*points):
    x = []
    y = []
    z = []
    for i in range(len(points)):
        x.append(points[i][0])
        y.append(points[i][1])
        z.append(points[i][2])
    return x,y,z

def line_plane_intersection(plane_p1, plane_p2, plane_p3, line_dir, line_point):
    plane_normal = np.cross(plane_p2 - plane_p1, plane_p3 - plane_p1)
    dot_prod = np.dot(plane_normal, line_dir)
    if abs(dot_prod) < 1e-8:
        return [0,0,0]
    v = line_point - plane_p1
    f = - np.dot(plane_normal, v) / dot_prod
    return v + f * line_dir + plane_p1

def perp_point(point_a1, point_b1, point_b2):
    point_b_norm = (point_b2 - point_b1) / np.linalg.norm(point_b2 - point_b1)
    v = point_a1 - point_b1
    t = np.dot(v, point_b_norm)
    return point_b1 + t * point_b_norm

def save_pdf(filename, label, size, segmentpos, segmenthoehe):
    
    if size == 'DIN A3':
        din_a4_width_mm  = 420
        din_a4_height_mm = 297
    elif size == 'DIN A5':
        din_a4_width_mm  = 210
        din_a4_height_mm = 148
    else: #DIN A4
        din_a4_width_mm  = 297
        din_a4_height_mm = 210
    
    din_a4_width_inches  = din_a4_width_mm  / 25.4
    din_a4_height_inches = din_a4_height_mm / 25.4
    ax_rel_width  = 0.9
    ax_rel_height = 0.9

    # figure with real units
    fig = plt.figure(figsize = (din_a4_width_inches, din_a4_height_inches))
    h = [Size.Fixed((1-ax_rel_width)  / 2 * din_a4_width_inches),  Size.Fixed(ax_rel_width  * din_a4_width_inches)]
    v = [Size.Fixed((1-ax_rel_height) / 2 * din_a4_height_inches), Size.Fixed(ax_rel_height * din_a4_height_inches)]
    divider = Divider(fig, (0, 0, 1, 1), h, v, aspect=False)
    ax = fig.add_axes(divider.get_position(),axes_locator=divider.new_locator(nx=1, ny=1))
    ax.set_xlim(-din_a4_width_mm  * ax_rel_width  / 2, din_a4_width_mm  * ax_rel_width  / 2)
    ax.set_ylim(-din_a4_height_mm * ax_rel_height * 0.05, din_a4_height_mm * ax_rel_height * 0.95)
    plt.gca().invert_yaxis()
    ax.xaxis.tick_top()

    # plot
    plt.plot(segmentpos, segmenthoehe, linewidth = 0.5, color = 'black') #, marker='.', markersize=0.01, markeredgecolor='black'
    plt.plot([0,0],[-din_a4_height_mm * ax_rel_height * 0.05, din_a4_height_mm * ax_rel_height * 0.95], linewidth=0.5, color='black', linestyle='--')
    plt.text(segmentpos[-2]+3,3,label,va='top',ha='left')

    # PDF speichern
    fig.savefig(filename)
    plt.close(fig)
    return

def calc_vns(mode, path, BG, L, D, H, Hmin, theta):

    ###############################################################
    # PARAMETER

    param = {}
    param['L']  = float(L)
    param['D']  = float(D)
    param['SW'] = np.arctan(param['D']/2/(param['L']))
    param['BG'] = np.radians(float(BG))
    param['H']  = float(H)
    param['Hmin'] = float(Hmin)
    param['J']  = param['L'] + param['D'] / 2 * np.tan(param['SW'])

    # Namen
    BASENAME = "Segment" +  "_BG-" + str(int(round(180 / np.pi * param['BG'] * 10, 0))) + \
                            "_L-" + str(int(param['L'])) + \
                            "_D-" + str(int(param['D'])) + \
                            "_H-" + str(int(param['H'])) + \
                            "_Hmin-" + str(int(param['Hmin']))
    print(BASENAME, "%.2f" % theta)

    ###############################################################
    # DEFINITION DER PLATTFORM

    vectors = {}
    vectors['suedlager'] = np.array([0, 0, param['H']])
    vectors['erdachse'] = vectors['suedlager'] + param['L'] * np.array([np.cos(param['BG']), 0, np.sin(param['BG'])])
    vectors['nordspitze']  = np.array([param['J'], 0,           param['H']])
    vectors['nordlot']     = np.array([param['J'], 0,           0])
    vectors['ostlager']    = np.array([param['L'], -param['D']/2, 0])
    vectors['westlager']   = np.array([param['L'], +param['D']/2, 0])
    vectors['ostsegment']  = np.array([param['L'], -param['D']/2, param['H']])
    vectors['westsegment'] = np.array([param['L'], +param['D']/2, param['H']])
    rot_vectors = copy.deepcopy(vectors)
    new_rot_vectors = copy.deepcopy(vectors)    
    
    
    
    ###############################################################
    # BERECHNUNG DER SEGMENTE

    segmentpos = []
    segmenthoehe = []
    segmentwinkel = []
    segmentgeschwindigkeit = []
    segmentlageroffset = []
    for th in np.arange(30, -30-0.3, -0.3):
        for key in vectors:
            if not key in ['suedlager','erdachse','ostlager','westlager']:
                rot_vectors[key] = vectors['suedlager'] + rot_around(vectors[key] - vectors['suedlager'], np.radians(th), vectors['erdachse'] - vectors['suedlager'])

        # Die folgenden zwei Zeilen ist der Kern der ganzen Berechnung
        # Es wird erst geschaut, wo der Lagerstrahl die Segmentebene schneidet, man bekommt einen Kontaktpunkt zwischen Lager und Segment
        # Anschließend wird vom Kontaktpunkt senkrecht auf die Segmentoberkante projiziert, um die Referenzposition des Kontakts für das Segment zu erhalten
        # Die Kontaktpunkte und Referenzpositionen die man bei verschiedenen Winkeln (th) erhält, definieren die Segmentform
        # Die Berechnung wird nur für die Ostseite durchgeführt, die Westseite ist gespiegelt gleich
        rot_vectors['ostkontakt'] = line_plane_intersection(rot_vectors['ostsegment'], rot_vectors['nordspitze'], rot_vectors['nordlot'], rot_vectors['ostlager'], rot_vectors['ostlager'])
        rot_vectors['ostref'] = perp_point(rot_vectors['ostkontakt'], rot_vectors['ostsegment'], rot_vectors['nordspitze'])
        this_segmentpos = np.linalg.norm(rot_vectors['ostref'] - rot_vectors['ostsegment'])
        this_segmenthoehe = np.linalg.norm(rot_vectors['ostkontakt'] - rot_vectors['ostref'])
        
        # Im Prinzip sind Referenzpositionen und Segmenthöhen einfache Abstände zwischen Punkten,
        # für eine stabile Berechnung auch in Extremfällen müssen jedoch die Vorzeichen ordentlich bestimmt werden
        if np.linalg.norm(rot_vectors['ostref'] - rot_vectors['nordspitze']) < np.linalg.norm(rot_vectors['ostsegment'] - rot_vectors['nordspitze']):
            this_segmentpos = -this_segmentpos
        if rot_vectors['ostref'][2] < 0:
            this_segmenthoehe = -this_segmenthoehe

        # Außen nicht über Plattformebene
        if not this_segmenthoehe > param['Hmin']:
            continue
        
        # Innen nicht über Symmetrieachse
        if abs(this_segmentpos) > np.linalg.norm(rot_vectors['ostsegment'] - rot_vectors['nordspitze']):
            continue

        # Innen nicht weiter als nötig (von außen vorgegeben)
        if segmentwinkel:
            if abs(th) > abs(segmentwinkel[0]):
                continue

        # Wenn nicht übersprungen, in Listen schreiben
        segmentpos.append(this_segmentpos)
        segmenthoehe.append(this_segmenthoehe)
        segmentwinkel.append(th)
        if not segmentgeschwindigkeit:
            segmentgeschwindigkeit.append(0)
        else:
            segmentgeschwindigkeit.append(np.sqrt((segmentpos[-2]-segmentpos[-1])**2 + (segmenthoehe[-2]-segmenthoehe[-1])**2) / (segmentwinkel[-2]-segmentwinkel[-1]) * 360 / 24)
        if np.linalg.norm(rot_vectors['ostkontakt']) > np.linalg.norm(rot_vectors['ostlager']):
            segmentlageroffset.append(np.linalg.norm(rot_vectors['ostkontakt'] - rot_vectors['ostlager']))
        else:
            segmentlageroffset.append(-np.linalg.norm(rot_vectors['ostkontakt'] - rot_vectors['ostlager']))

        # Erster geschriebener Punkt definiert die Ecken
        if len(segmentpos) == 1:
            new_rot_vectors['osteck'] = rot_vectors['ostref']
            new_rot_vectors['westeck'] = rot_vectors['nordspitze'] + np.linalg.norm(new_rot_vectors['osteck'] - rot_vectors['nordspitze']) * (rot_vectors['westsegment'] - rot_vectors['nordspitze']) / np.linalg.norm(rot_vectors['westsegment'] - rot_vectors['nordspitze'])
        
        # Letzter geschriebener Punkt definiert die Innenkanten (einfach jedesmal überschreiben)
        new_rot_vectors['ostinnen'] = rot_vectors['ostref']
        new_rot_vectors['westinnen'] = rot_vectors['nordspitze'] + np.linalg.norm(new_rot_vectors['ostinnen'] - rot_vectors['nordspitze']) * (rot_vectors['westsegment'] - rot_vectors['nordspitze']) / np.linalg.norm(rot_vectors['westsegment'] - rot_vectors['nordspitze'])
    
    if not segmentwinkel:
        print("Kein einziger Punkt des Segments wurde gezeichnet!")
        print("Interrupt!")
        return

    # Rotiere gemerkte Punkte zurück
    vectors['osteck']    = rot_vectors['suedlager'] + rot_around(new_rot_vectors['osteck']    - rot_vectors['suedlager'], -np.radians(segmentwinkel[0]),  rot_vectors['erdachse'] - rot_vectors['suedlager'])
    vectors['westeck']   = rot_vectors['suedlager'] + rot_around(new_rot_vectors['westeck']   - rot_vectors['suedlager'], -np.radians(segmentwinkel[0]),  rot_vectors['erdachse'] - rot_vectors['suedlager'])
    vectors['ostinnen']  = rot_vectors['suedlager'] + rot_around(new_rot_vectors['ostinnen']  - rot_vectors['suedlager'], -np.radians(segmentwinkel[-1]), rot_vectors['erdachse'] - rot_vectors['suedlager'])
    vectors['westinnen'] = rot_vectors['suedlager'] + rot_around(new_rot_vectors['westinnen'] - rot_vectors['suedlager'], -np.radians(segmentwinkel[-1]), rot_vectors['erdachse'] - rot_vectors['suedlager'])

    # Anfang und Ende der Segmente auf Höhe Null
    segmentpos.insert(0, segmentpos[0])
    segmentpos.append(segmentpos[-1])
    segmentpos.append(segmentpos[0])
    segmenthoehe.insert(0, 0)
    segmenthoehe.append(0)
    segmenthoehe.append(0)
    segmentgeschwindigkeit[0] = 2 * segmentgeschwindigkeit[1] - segmentgeschwindigkeit[2]
    segmentgeschwindigkeit_norm = [x/min(segmentgeschwindigkeit) for x in segmentgeschwindigkeit]

    # Nur return
    if mode==0:
        return BASENAME, segmentwinkel[0]

    # Ausgeben und abspeichern
    if mode==3 and not not path:
        with open(os.path.join(path, BASENAME + "_Segment.dat"), 'w') as f:
            for i in range(len(segmentpos)):
                f.write("%.6f" % segmentpos[i] + "\t")
                f.write("%.6f" % -segmenthoehe[i] + "\n")
        with open(os.path.join(path, BASENAME + "_Geschwindigkeit.dat"), 'w') as f:
            for i in range(len(segmentwinkel)):
                f.write("%.6f" % segmentwinkel[i] + "\t")
                f.write("%.6f" % segmentgeschwindigkeit[i] + "\n")
        with open(os.path.join(path, BASENAME + "_Geschwindigkeit_norm.dat"), 'w') as f:
            for i in range(len(segmentwinkel)):
                f.write("%.6f" % segmentwinkel[i] + "\t")
                f.write("%.6f" % segmentgeschwindigkeit_norm[i] + "\n")
        with open(os.path.join(path, BASENAME + "_Lageroffset.dat"), 'w') as f:
            for i in range(len(segmentwinkel)):
                f.write("%.6f" % segmentwinkel[i] + "\t")
                f.write("%.6f" % segmentlageroffset[i] + "\n")

        mpl.use("pdf")
        save_pdf(os.path.join(path, BASENAME  + "_DIN_A3.pdf"), BASENAME.replace("_",",  ").replace("-",": ").replace("Segment,  ",""), 'DIN A3', segmentpos, segmenthoehe)
        save_pdf(os.path.join(path, BASENAME  + "_DIN_A4.pdf"), BASENAME.replace("_",",  ").replace("-",": ").replace("Segment,  ",""), 'DIN A4', segmentpos, segmenthoehe)
        save_pdf(os.path.join(path, BASENAME  + "_DIN_A5.pdf"), BASENAME.replace("_",",  ").replace("-",": ").replace("Segment,  ",""), 'DIN A5', segmentpos, segmenthoehe)


    ###############################################################
    # Plattform Ansicht

    if mode in [1,2]:
        
        if mode == 2:
            mpl.use("Agg")
        else:
            mpl.use("Qt5Agg")
        
        theta = theta + 0.00001 # vertical surfaces impossible to plot

        # Gehe zu Zeichnen-Position
        for key in vectors:
            if not key in ['suedlager','erdachse','ostlager','westlager']:
                rot_vectors[key] = vectors['suedlager'] + rot_around(vectors[key] - vectors['suedlager'], np.radians(theta), vectors['erdachse'] - vectors['suedlager'])
        rot_vectors['ostkontakt'] = line_plane_intersection(rot_vectors['ostsegment'], rot_vectors['nordspitze'], rot_vectors['nordlot'], rot_vectors['ostlager'], rot_vectors['ostlager'])
        rot_vectors['westkontakt'] = line_plane_intersection(rot_vectors['westsegment'], rot_vectors['nordspitze'], rot_vectors['nordlot'], rot_vectors['westlager'], rot_vectors['westlager'])
        rot_vectors['ostref'] = perp_point(rot_vectors['ostkontakt'], rot_vectors['osteck'], rot_vectors['nordspitze'])
        rot_vectors['westref'] = perp_point(rot_vectors['westkontakt'], rot_vectors['westeck'], rot_vectors['nordspitze'])    

        # Erzeuge Figure
        fig, ax = plt.subplots(subplot_kw={"projection": "3d"}, figsize = (5, 5))

        # Plattformebene
        [x,y,z] = points_to_coordinates(
            rot_vectors['suedlager'],
            rot_vectors['osteck'],
            rot_vectors['ostinnen'],
            rot_vectors['westinnen'],
            rot_vectors['westeck'],
            )
        ax.plot_trisurf(x,y,z, color=[0.6,0.6,0.6], alpha=0.6, shade=False, antialiased=False)

        # Plattform-Kanten
        [x,y,z] = points_to_coordinates([0,0,0], rot_vectors['suedlager'])
        ax.plot(x,y,z, linewidth=1, color=[0.2,0.2,0.2])
        [x,y,z] = points_to_coordinates(rot_vectors['suedlager'], rot_vectors['osteck'])
        ax.plot(x,y,z, linewidth=1, color=[0.2,0.2,0.2])
        [x,y,z] = points_to_coordinates(rot_vectors['suedlager'], rot_vectors['westeck'])
        ax.plot(x,y,z, linewidth=1, color=[0.2,0.2,0.2])
        [x,y,z] = points_to_coordinates(rot_vectors['osteck'], rot_vectors['ostinnen'])
        ax.plot(x,y,z, linewidth=1, color=[0.2,0.2,0.2])
        [x,y,z] = points_to_coordinates(rot_vectors['westeck'], rot_vectors['westinnen'])
        ax.plot(x,y,z, linewidth=1, color=[0.2,0.2,0.2])
        [x,y,z] = points_to_coordinates(rot_vectors['ostinnen'], rot_vectors['westinnen'])
        ax.plot(x,y,z, linewidth=1, color=[0.2,0.2,0.2])

        # Erdachse
        [x,y,z] = points_to_coordinates(rot_vectors['suedlager'], rot_vectors['erdachse'])
        ax.plot(x,y,z, linewidth=1, color='red')

        # Nordlager
        [x,y,z] = points_to_coordinates([0,0,0], 1.1 * rot_vectors['ostlager'])
        ax.plot(x,y,z, linewidth=1, color='blue')
        [x,y,z] = points_to_coordinates([0,0,0], 1.1 * rot_vectors['westlager'])
        ax.plot(x,y,z, linewidth=1, color='blue')

        # Segmentebenen
        if False:
            [x,y,z] = points_to_coordinates(
                rot_vectors['nordspitze'],
                rot_vectors['osteck'] + 0.2 * (rot_vectors['osteck'] - rot_vectors['nordspitze']),
                rot_vectors['osteck'] + 0.2 * (rot_vectors['osteck'] - rot_vectors['nordspitze']) + 2 * (rot_vectors['nordlot'] - rot_vectors['nordspitze']),
                rot_vectors['nordspitze'] + 2 * (rot_vectors['nordlot'] - rot_vectors['nordspitze']),
                )
            ax.plot_trisurf(x,y,z, color='green', alpha=0.4)
            [x,y,z] = points_to_coordinates(
                rot_vectors['nordspitze'],
                rot_vectors['westeck'] + 0.2 * (rot_vectors['westeck'] - rot_vectors['nordspitze']),
                rot_vectors['westeck'] + 0.2 * (rot_vectors['westeck'] - rot_vectors['nordspitze']) + 2 * (rot_vectors['nordlot'] - rot_vectors['nordspitze']),
                rot_vectors['nordspitze'] + 2 * (rot_vectors['nordlot'] - rot_vectors['nordspitze']),
                )
            ax.plot_trisurf(x,y,z, color='green', alpha=0.4)
            [x,y,z] = points_to_coordinates(rot_vectors['ostkontakt'], rot_vectors['ostref'])
            ax.plot(x,y,z, linewidth=1, color='green')
            [x,y,z] = points_to_coordinates(rot_vectors['westkontakt'], rot_vectors['westref'])
            ax.plot(x,y,z, linewidth=1, color='green')

        # Segmente
        else:

            # Ostsegment
            segmentpunkte_x = []
            segmentpunkte_y = []
            segmentpunkte_z = []
            for i in range(len(segmentpos)):
                ref = rot_vectors['ostsegment'] + segmentpos[i] * (rot_vectors['ostsegment'] - rot_vectors['nordspitze']) / np.linalg.norm(rot_vectors['ostsegment'] - rot_vectors['nordspitze'])
                punkt = ref + segmenthoehe[i] * (rot_vectors['nordlot'] - rot_vectors['nordspitze']) / np.linalg.norm(rot_vectors['nordlot'] - rot_vectors['nordspitze'])
                segmentpunkte_x.append(punkt[0])
                segmentpunkte_y.append(punkt[1])
                segmentpunkte_z.append(punkt[2])
            ax.plot_trisurf(segmentpunkte_x,segmentpunkte_y,segmentpunkte_z, color=[0.2,0.6,0.2], alpha=1, shade=False, antialiased=False)
            ax.plot(segmentpunkte_x,segmentpunkte_y,segmentpunkte_z, linewidth=1, color='green')

            # Westsegment
            segmentpunkte_x = []
            segmentpunkte_y = []
            segmentpunkte_z = []
            for i in range(len(segmentpos)):
                ref = rot_vectors['westsegment'] + segmentpos[i] * (rot_vectors['westsegment'] - rot_vectors['nordspitze']) / np.linalg.norm(rot_vectors['westsegment'] - rot_vectors['nordspitze'])
                punkt = ref + segmenthoehe[i] * (rot_vectors['nordlot'] - rot_vectors['nordspitze']) / np.linalg.norm(rot_vectors['nordlot'] - rot_vectors['nordspitze'])
                segmentpunkte_x.append(punkt[0])
                segmentpunkte_y.append(punkt[1])
                segmentpunkte_z.append(punkt[2])
            ax.plot_trisurf(segmentpunkte_x,segmentpunkte_y,segmentpunkte_z, color=[0.2,0.6,0.2], alpha=1, shade=False, antialiased=False)
            ax.plot(segmentpunkte_x,segmentpunkte_y,segmentpunkte_z, linewidth=1, color='green')

        # Limits
        boxwidth = param['L']*1.2
        ax.set_xlim(-5, -5+boxwidth)
        ax.set_ylim(-boxwidth/2, boxwidth/2)
        ax.set_zlim(0, boxwidth)

        # Lagerebene
        [x,y,z] = points_to_coordinates(    
            [ax.get_xlim()[0],ax.get_ylim()[0],0],
            [ax.get_xlim()[0],ax.get_ylim()[1],0],
            [ax.get_xlim()[1],ax.get_ylim()[0],0],
            [ax.get_xlim()[1],ax.get_ylim()[1],0]
            )
        ax.plot_trisurf(x,y,z, color='grey', alpha=0.1, shade=False)

        # Anzeige
        ax.set_axis_off()
        ax.view_init(22,-30)

    if mode == 1:
        plt.show()

    if mode == 2:
        fig.canvas.draw()       # draw the canvas, cache the renderer
        image = np.frombuffer(fig.canvas.tostring_rgb(), dtype='uint8')
        print("reshape " + str(len(image)) + " to " + str(fig.canvas.get_width_height()[::-1] + (3,)))
        image = image.reshape(fig.canvas.get_width_height()[::-1] + (3,))
        plt.close()
        return image


class NewGUI():
    def __init__(self):
        
        self.root = tk.Tk()
        self.root.title("VNS-Plattform (" + GUI_version + ")")
        self.root.resizable(width=False, height=False)
        self.complete = 0

        # reset saved window position and path
        if os.path.isfile("VNS_Plattform.conf"): 
            with open("VNS_Plattform.conf", "r") as conf:
                lines = conf.readlines()
                self.restore = lines[0].split('#')
        else:
            self.restore = []
        self.root.protocol("WM_DELETE_WINDOW",  self.on_close)

        try:
            self.root.iconbitmap('VNS.ico')
        except:
            print("Icon nicht gefunden... Fahre fort.")

        if not os.path.isfile("VNS_Plattform_Skizze.png"):
            print("Skizze nicht gefunden!! ... Fahre fort")
        else:
            pix_per_pts = self.root.winfo_fpixels('1p')
            print("%.2f" % pix_per_pts + " pixels per point")
            img = ImageTk.PhotoImage(Image.open('VNS_Plattform_Skizze.png').resize((int(160*pix_per_pts),int(240*pix_per_pts)), Image.LANCZOS))
            image = tk.Label(image=img)
            image.grid(row=0, column=0, rowspan=8, columnspan=2)

        padx = 10
        pady = 5
    
        BG_row = 0
        self.label_BG = tk.Label(text="Breitengrad BG (°)")
        self.label_BG.grid(row=BG_row, column=2, sticky='w')
        self.entry_BG = tk.Entry()
        self.entry_BG.grid(row=BG_row, column=3, padx=padx, pady=pady, sticky='news')
        self.result_BG = tk.Entry(state='disabled')
        self.result_BG.grid(row=BG_row, column=4, padx=padx, pady=pady, sticky='news')

        S_row = 1
        self.label_S = tk.Label(text="Schwerpunkt S (mm)")
        self.label_S.grid(row=S_row, column=2, sticky='w')
        self.entry_S = tk.Entry()
        self.entry_S.grid(row=S_row, column=3, padx=padx, pady=pady, sticky='news')
        self.result_S = tk.Entry(state='disabled')
        self.result_S.grid(row=S_row, column=4, padx=padx, pady=pady, sticky='news')

        B_row = 2
        self.label_L = tk.Label(text="Länge L (mm)")
        self.label_L.grid(row=B_row, column=2, sticky='w')
        self.entry_L = tk.Entry()
        self.entry_L.grid(row=B_row, column=3, padx=padx, pady=pady, sticky='news')
        self.result_L = tk.Entry(state='disabled')
        self.result_L.grid(row=B_row, column=4, padx=padx, pady=pady, sticky='news')

        F_row = 3
        self.label_D = tk.Label(text="Breite D (mm)")
        self.label_D.grid(row=F_row, column=2, sticky='w')
        self.entry_D = tk.Entry()
        self.entry_D.grid(row=F_row, column=3, padx=padx, pady=pady, sticky='news')
        self.result_D = tk.Entry(state='disabled')
        self.result_D.grid(row=F_row, column=4, padx=padx, pady=pady, sticky='news')        

        SW_row = 4
        self.label_SW = tk.Label(text="Segmentwinkel SW (°)")
        self.label_SW.grid(row=SW_row, column=2, sticky='w')
        self.entry_SW = tk.Entry()
        self.entry_SW.grid(row=SW_row, column=3, padx=padx, pady=pady, sticky='news')
        self.result_SW = tk.Entry(state='disabled')
        self.result_SW.grid(row=SW_row, column=4, padx=padx, pady=pady, sticky='news')

        H_row = 5
        self.label_H = tk.Label(text="Höhe H (mm)")
        self.label_H.grid(row=H_row, column=2, sticky='w')
        self.entry_H = tk.Entry()
        self.entry_H.grid(row=H_row, column=3, padx=padx, pady=pady, sticky='news')
        self.result_H = tk.Entry(state='disabled')
        self.result_H.grid(row=H_row, column=4, padx=padx, pady=pady, sticky='news')

        Hmin_row = 6
        self.label_Hmin = tk.Label(text="Mindesthöhe Hmin (mm)")
        self.label_Hmin.grid(row=Hmin_row, column=2, sticky='w')
        self.entry_Hmin = tk.Entry()
        self.entry_Hmin.grid(row=Hmin_row, column=3, padx=padx, pady=pady, sticky='news')
        self.result_Hmin = tk.Entry(state='disabled')
        self.result_Hmin.grid(row=Hmin_row, column=4, padx=padx, pady=pady, sticky='news')
        
        button_row = 7
        self.button_show = tk.Button(text="Anzeigen", command=self.show)
        self.button_show.grid(row=button_row, column=2, padx=padx, pady=pady, sticky="news")
        self.button_gif = tk.Button(text="GIF", command=self.gif)
        self.button_gif.grid(row=button_row, column=3, padx=padx, pady=pady, sticky="news")
        self.button_process = tk.Button(text="PDF & DATA", command=self.save)
        self.button_process.grid(row=button_row, column=4, padx=padx, pady=pady, sticky="news")

        self.entries = [self.entry_BG, self.entry_S, self.entry_L, self.entry_D, self.entry_SW, self.entry_H, self.entry_Hmin]
        self.results = [self.result_BG, self.result_S, self.result_L, self.result_D, self.result_SW, self.result_H, self.result_Hmin]

        for i in range(len(self.entries)):
            self.entries[i].bind('<KeyRelease>', self.check_entries)
            self.entries[i].bind('<FocusOut>', self.check_entries)
            if len(self.restore)==7:
                set_entry(self.entries[i], self.restore[i])

        self.check_entries()
        self.root.mainloop()
        
    def show(self):
        try:
            if self.complete:
                calc_vns(1, '', self.result_BG.get(), self.result_L.get(), self.result_D.get(), self.result_H.get(), self.result_Hmin.get(), 0)
        except Exception as e:
            print(e)
            tk.messagebox.showerror('Error', 'Unknown error!')
            return

    def save(self):
        try:
            path = tk.filedialog.askdirectory(title="Choose path")
            print(path)
            if path and self.complete:
                calc_vns(3, path, self.result_BG.get(), self.result_L.get(), self.result_D.get(), self.result_H.get(), self.result_Hmin.get(), 0)
        except Exception as e:
            print(e)
            tk.messagebox.showerror('Error', 'Unknown error!')
            return

    def gif(self):
        try:
            vns_info = calc_vns(0, '', self.result_BG.get(), self.result_L.get(), self.result_D.get(), self.result_H.get(), self.result_Hmin.get(), 0)
            filename = tk.filedialog.asksaveasfilename(initialfile=os.path.join(vns_info[0] + '.gif'), title="Choose filename", defaultextension=".gif", filetypes=(("Graphics Interchange Format", "*.gif"),("All Files", "*.*")))
            print(filename)
            if filename and self.complete:
                print("imageio version " + imageio.__version__)
                imageio.mimsave(filename, [calc_vns(2, '', self.result_BG.get(), self.result_L.get(), self.result_D.get(), self.result_H.get(), self.result_Hmin.get(), th)
                                 for th in np.arange(-vns_info[1], vns_info[1], 2*vns_info[1]/20)], fps=10)
        except Exception as e:
            print(e)
            tk.messagebox.showerror('Error', 'Unknown error!')
            return

    def check_entries(self, *args):
        for i in range(len(self.entries)):
            process_entry(self.entries[i], self.results[i])
        
        BG   = self.result_BG.get()
        S    = self.result_S.get()
        B    = self.result_L.get()
        F    = self.result_D.get()
        SW   = self.result_SW.get()
        H    = self.result_H.get()
        Hmin = self.result_Hmin.get()  
        
        # check dependencies
        if not BG:
            set_result(self.result_BG, "Def. <BG>")
            set_result(self.result_S,  "Def. <BG>")
            set_result(self.result_L,  "Def. <BG>")
        else:
            if not S:
                if not B:
                    set_result(self.result_S, "Def. <S> / <L>")
                    set_result(self.result_L, "Def. <S> / <L>")
                else:
                    S = "%.2f" % (float(B) * np.sin(np.radians(float(BG))) * np.cos(np.radians(float(BG))))
                    set_result(self.result_S, S)
            elif not B:
                B = "%.2f" % (float(S) / np.sin(np.radians(float(BG))) / np.cos(np.radians(float(BG))))
                set_result(self.result_L, B)
            else:
                set_result(self.result_S, "!! <S> / <L> !!")
                set_result(self.result_L, "!! <S> / <L> !!")
        if not B:
            set_result(self.result_D,  "Def. <L>")
            set_result(self.result_SW, "Def. <L>")
        else:
            if not F:
                if not SW:
                    set_result(self.result_D,  "Def. <D> / <SW>")
                    set_result(self.result_SW, "Def. <D> / <SW>")
                else:
                    F = "%.2f" % (2 * float(B) * np.tan(np.radians(float(SW))))
                    set_result(self.result_D, F)
            elif not SW:
                SW = "%.2f" % (np.degrees(np.arctan(float(F)/2/float(B))))
                set_result(self.result_SW, SW)
            else:
                set_result(self.result_D,  "!! <D> / <SW> !!")
                set_result(self.result_SW, "!! <D> / <SW> !!")
        if not H:
            set_result(self.result_H, "Def. <H>")
        if not Hmin:
            set_result(self.result_Hmin, "Def. <Hmin>")
            
        # check complete
        if BG and S and B and F and SW and H and Hmin:
            self.complete = 1
        else:
            self.complete = 0

    def on_close(self):
        with open("VNS_Plattform.conf", "w") as conf: 
            conf.write(    self.entry_BG.get() + "#"
                         + self.entry_S.get() + "#"
                         + self.entry_L.get() + "#"
                         + self.entry_D.get() + "#"
                         + self.entry_SW.get() + "#"
                         + self.entry_H.get() + "#"
                         + self.entry_Hmin.get())
        self.root.destroy()


if __name__ == '__main__':
    new = NewGUI()

