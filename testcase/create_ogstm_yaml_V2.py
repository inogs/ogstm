import os
import pyfabm
import tkinter as tk
from tkinter import ttk, messagebox

CODEPATH = '../../'
CODEPATH = CODEPATH.replace("~", os.getenv("HOME"))
#fabm_yaml = CODEPATH + "/fabm/extern/ogs/fabm_multispectral_2xDetritus.yaml"
fabm_yaml = "fabm_rosenmcartur.yaml"
class OGSTMYAMLGenerator:
    def __init__(self, root):
        self.root = root
        self.root.title("OGSTM YAML Generator")
        self.root.geometry("800x600")
        
        # Load FABM model
        try:
            self.model = pyfabm.Model(fabm_yaml)
        except Exception as e:
            messagebox.showerror("Error", f"Failed to load FABM model: {e}")
            return
        
        self.create_widgets()
        
    def create_widgets(self):
        # Main frame with scrollbar
        main_frame = ttk.Frame(self.root)
        main_frame.pack(fill=tk.BOTH, expand=True, padx=10, pady=10)
        
        # Canvas and scrollbar for scrollable content
        canvas = tk.Canvas(main_frame)
        scrollbar = ttk.Scrollbar(main_frame, orient="vertical", command=canvas.yview)
        scrollable_frame = ttk.Frame(canvas)
        
        scrollable_frame.bind(
            "<Configure>",
            lambda e: canvas.configure(scrollregion=canvas.bbox("all"))
        )
        
        canvas.create_window((0, 0), window=scrollable_frame, anchor="nw")
        canvas.configure(yscrollcommand=scrollbar.set)
        
        # Interior Diagnostics Section
        interior_label = ttk.Label(scrollable_frame, text="Interior Diagnostic Variables", 
                                 font=("Arial", 12, "bold"))
        interior_label.pack(anchor="w", pady=(0, 10))
        
        self.interior_vars = {}
        interior_frame = ttk.Frame(scrollable_frame)
        interior_frame.pack(fill="x", padx=10, pady=(0, 20))
        
        # Headers
        ttk.Label(interior_frame, text="Variable Name", font=("Arial", 10, "bold")).grid(
            row=0, column=0, sticky="w", padx=(0, 20))
        ttk.Label(interior_frame, text="diahf", font=("Arial", 10, "bold")).grid(
            row=0, column=1, padx=10)
        ttk.Label(interior_frame, text="diaWR", font=("Arial", 10, "bold")).grid(
            row=0, column=2, padx=10)
        
        row = 1
        for variable in self.model.interior_diagnostic_variables:
            if variable.output:
                var_name = variable.name.replace('/', '_')
                
                # Variable name
                ttk.Label(interior_frame, text=var_name).grid(
                    row=row, column=0, sticky="w", padx=(0, 20), pady=2)
                
                # diahf checkbox
                diahf_var = tk.BooleanVar()
                ttk.Checkbutton(interior_frame, variable=diahf_var).grid(
                    row=row, column=1, padx=10, pady=2)
                
                # diaWR checkbox
                diawr_var = tk.BooleanVar()
                ttk.Checkbutton(interior_frame, variable=diawr_var).grid(
                    row=row, column=2, padx=10, pady=2)
                
                self.interior_vars[var_name] = {'diahf': diahf_var, 'diaWR': diawr_var}
                row += 1
        
        # Horizontal Diagnostics Section
        horizontal_label = ttk.Label(scrollable_frame, text="Horizontal Diagnostic Variables", 
                                   font=("Arial", 12, "bold"))
        horizontal_label.pack(anchor="w", pady=(20, 10))
        
        self.horizontal_vars = {}
        horizontal_frame = ttk.Frame(scrollable_frame)
        horizontal_frame.pack(fill="x", padx=10, pady=(0, 20))
        
        # Headers
        ttk.Label(horizontal_frame, text="Variable Name", font=("Arial", 10, "bold")).grid(
            row=0, column=0, sticky="w", padx=(0, 20))
        ttk.Label(horizontal_frame, text="diahf_2d", font=("Arial", 10, "bold")).grid(
            row=0, column=1, padx=10)
        ttk.Label(horizontal_frame, text="diaWR_2d", font=("Arial", 10, "bold")).grid(
            row=0, column=2, padx=10)
        
        row = 1
        for variable in self.model.horizontal_diagnostic_variables:
            if variable.output:
                var_name = variable.name.replace('/', '_')
                
                # Variable name
                ttk.Label(horizontal_frame, text=var_name).grid(
                    row=row, column=0, sticky="w", padx=(0, 20), pady=2)
                
                # diahf_2d checkbox
                diahf_2d_var = tk.BooleanVar()
                ttk.Checkbutton(horizontal_frame, variable=diahf_2d_var).grid(
                    row=row, column=1, padx=10, pady=2)
                
                # diaWR_2d checkbox
                diawr_2d_var = tk.BooleanVar()
                ttk.Checkbutton(horizontal_frame, variable=diawr_2d_var).grid(
                    row=row, column=2, padx=10, pady=2)
                
                self.horizontal_vars[var_name] = {'diahf_2d': diahf_2d_var, 'diaWR_2d': diawr_2d_var}
                row += 1
        
        # Pack canvas and scrollbar
        canvas.pack(side="left", fill="both", expand=True)
        scrollbar.pack(side="right", fill="y")
        
        # Buttons frame
        button_frame = ttk.Frame(self.root)
        button_frame.pack(fill="x", padx=10, pady=10)
        
        ttk.Button(button_frame, text="Generate YAML", 
                  command=self.generate_yaml).pack(side="right", padx=(10, 0))
        ttk.Button(button_frame, text="Select All", 
                  command=self.select_all).pack(side="right", padx=(10, 0))
        ttk.Button(button_frame, text="Deselect All", 
                  command=self.deselect_all).pack(side="right")
        
        # Bind mousewheel to canvas
        def on_mousewheel(event):
            canvas.yview_scroll(int(-1*(event.delta/120)), "units")
        canvas.bind_all("<MouseWheel>", on_mousewheel)
    
    def select_all(self):
        for var_dict in self.interior_vars.values():
            var_dict['diahf'].set(True)
            var_dict['diaWR'].set(True)
        for var_dict in self.horizontal_vars.values():
            var_dict['diahf_2d'].set(True)
            var_dict['diaWR_2d'].set(True)
    
    def deselect_all(self):
        for var_dict in self.interior_vars.values():
            var_dict['diahf'].set(False)
            var_dict['diaWR'].set(False)
        for var_dict in self.horizontal_vars.values():
            var_dict['diahf_2d'].set(False)
            var_dict['diaWR_2d'].set(False)
    
    def generate_yaml(self):
        try:
            with open('ogstm.yaml', 'w') as f:
                # Interior state variables
                f.write('interior_state:\n')
                for i, variable in enumerate(self.model.state_variables):
                    f.write(f"    {variable.name.replace('/','_')}:\n")
                    f.write(f"        ctrmax: 1.000000e+03\n")
                    # OGSTM want always at least an high freq output
                    if i == 0:
                        f.write(f"        ctrhf: 1\n")
                    else:
                        f.write(f"        ctrhf: 0\n")
                    f.write(f"        relax: 0\n")
                
                # Interior diagnostic variables - only include if diaWR is True
                f.write('interior_diagnostic:\n')
                interior_written = False
                for var_name, var_dict in self.interior_vars.items():
                    if var_dict['diaWR'].get():  # Only write if diaWR is checked
                        f.write(f"    {var_name}:\n")
                        f.write(f"        diahf: {1 if var_dict['diahf'].get() else 0}\n")
                        f.write(f"        diaWR: 1\n")
                        interior_written = True
                
                # If no interior diagnostic variables were written, add a placeholder comment
                if not interior_written:
                    f.write('    # No interior diagnostic variables selected for output\n')
                
                # Horizontal diagnostic variables - only include if diaWR_2d is True
                f.write('horizontal_diagnostic:\n')
                horizontal_written = False
                for var_name, var_dict in self.horizontal_vars.items():
                    if var_dict['diaWR_2d'].get():  # Only write if diaWR_2d is checked
                        f.write(f"    {var_name}:\n")
                        f.write(f"        diahf_2d: {1 if var_dict['diahf_2d'].get() else 0}\n")
                        f.write(f"        diaWR_2d: 1\n")
                        horizontal_written = True
                
                # If no horizontal diagnostic variables were written, add a placeholder comment
                if not horizontal_written:
                    f.write('    # No horizontal diagnostic variables selected for output\n')
            
            messagebox.showinfo("Success", "ogstm.yaml file generated successfully!")
            
        except Exception as e:
            messagebox.showerror("Error", f"Failed to generate YAML file: {e}")

if __name__ == "__main__":
    root = tk.Tk()
    app = OGSTMYAMLGenerator(root)
    root.mainloop()
