import tkinter as tk
from tkinter import ttk, filedialog

FIELDS = [
    # (label, key, default, required, is_file)
    ("--- Required ---",      None,                    None,                                        None,  None), # label
    ("Input File 1 *",        "infile1",               "",                                          True,  True),
    ("Input File 2 *",        "infile2",               "",                                          True,  True),
    ("ID Column *",           "id_col",                "analyte_id",                                True,  False),
    ("Numeric Column *",      "numeric_col",           "Abundance rate",                            True,  False),
    ("--- Optional ---",      None,                    None,                                        None,  None), # label
    ("Only Common Genes",     "only_common_genes",     True,                                        False, None),  # bool
    ("Top N Annotations",     "n_top_annotations",     10,                                          False, False),
    ("Knowledgebase Dir",     "knowledgebase_directory","knowledgebase",                            False, True),
    ("Annotation Outfile 1",  "annotation_outfile1",   "outputs/functional_annotation_results1.tsv",False, False),
    ("Annotation Outfile 2",  "annotation_outfile2",   "outputs/functional_annotation_results_2.tsv",False,False),
    ("Shared Outfile",        "shared_outfile",        "outputs/shared_annotation_genes.tsv",       False, False),
    ("Final Outfile",         "final_outfile",         "outputs/final_comparison_results.tsv",      False, False),
]

def run_gui():
    result = {}
    root = tk.Tk()
    root.title("Functional Annotation Settings")
    root.resizable(False, False)

    frame = ttk.Frame(root, padding=10)
    frame.grid()

    vars_ = {}
    for row, (label, key, default, required, is_file) in enumerate(FIELDS):
        if key is None:  # separator
            ttk.Separator(frame, orient="horizontal").grid(row=row, column=0, columnspan=3, sticky="ew", pady=6)
            ttk.Label(frame, text=label, foreground="gray").grid(row=row, column=0, columnspan=3)
            continue
        ttk.Label(frame, text=label).grid(row=row, column=0, sticky="w", pady=2)

        if isinstance(default, bool):  # checkbox
            v = tk.BooleanVar(value=default)
            ttk.Checkbutton(frame, variable=v).grid(row=row, column=1, sticky="w")
        else:
            v = tk.StringVar(value=str(default))
            ttk.Entry(frame, textvariable=v, width=45).grid(row=row, column=1, pady=2)
            if is_file:
                cmd = (lambda v=v: v.set(filedialog.askopenfilename())) if is_file == True and "Dir" not in label \
                      else (lambda v=v: v.set(filedialog.askdirectory()))
                ttk.Button(frame, text="Browse", command=cmd).grid(row=row, column=2, padx=4)

        vars_[key] = (v, type(default))

    def submit():
        for label, key, default, required, _ in FIELDS:
            if key is None:
                continue
            v, typ = vars_[key]
            val = v.get()
            if required and not val.strip():
                tk.messagebox.showerror("Missing", f"{label.strip('* ')} is required.")
                return
            result[key] = typ(val) if typ != bool else v.get()
        root.destroy()

    ttk.Button(frame, text="Run", command=submit).grid(row=len(FIELDS), column=1, pady=10)
    root.mainloop()
    return result if result else None


if __name__ == "__main__":
    params = run_gui()
    if params:
        print(params)
