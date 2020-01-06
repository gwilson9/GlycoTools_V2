using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Data;
using System.Drawing;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows.Forms;
using System.Data.SQLite;
using System.IO;
using System.Reflection;
using System.Threading;
using ZedGraph;
using CSMSL;
using CSMSL.Chemistry;
using CSMSL.Proteomics;
using CSMSL.Spectral;
using CSMSL.IO.Thermo;
using LiveCharts;
using LiveCharts.Wpf;
using LiveCharts.WinForms;
using LiveCharts.Defaults;
using LiveCharts.Events;
using System.Globalization;
using System.Windows.Media;
using System.Windows.Media.Imaging;
using System.Windows.Controls;
using System.Diagnostics;

namespace _20190618_GlycoTools_V2
{
    public partial class GlycoTools : Form
    {
        //Data tables
        private BindingSource bindingSource1 = new BindingSource();
        //View peptides table
        private BindingSource bindingSource2 = new BindingSource();
        private SQLiteDataAdapter dataAdapter = new SQLiteDataAdapter();
        private bool dataUploaded = false;
        private BindingSource inSourceFragBindingSource = new BindingSource();
        private BindingSource dataUploadBindingSouce = new BindingSource();
        private bool dataFromDB = false;

        public GlycoTools()
        {
            InitializeComponent();
            this.DragEnter += new DragEventHandler(Form1_DragEnter);
            this.DragDrop += new DragEventHandler(Form1_DragDrop);
            var comboboxOptions = new int[] { 1, 2, 3 };
            this.modCountFilter.DataSource = comboboxOptions;
            this.SetStyle(ControlStyles.ResizeRedraw, true);
            dataTree.Tag = false;
            InitializePlots();
        }

        private void InitializePlots()
        {
            inSourceFragData.Dock = DockStyle.Fill;
            inSourceFragData.DataSource = inSourceFragBindingSource;
            inSourceFragChart.Dock = DockStyle.Fill;

            splitContainer7.Dock = DockStyle.Fill;

            inSourceFragChart.AxisY.Clear();
            inSourceFragChart.AxisY.Add(new LiveCharts.Wpf.Axis
            {
                Title = "Intensity",
                MinValue = 0.0,
                FontSize = 15,
                Foreground = System.Windows.Media.Brushes.Black,
                Sections = new SectionsCollection
                {
                    new AxisSection
                    {
                        Value = 0.00001,
                        StrokeThickness = 1,
                        Stroke = new SolidColorBrush(System.Windows.Media.Color.FromRgb(64,79,86))
                    }
                }
            });

            inSourceFragChart.AxisX.Clear();
            inSourceFragChart.AxisX.Add(new LiveCharts.Wpf.Axis
            {
                Title = "Time (min.)",
                FontSize = 15,
                Foreground = System.Windows.Media.Brushes.Black,
                Sections = new SectionsCollection
                {
                    new AxisSection
                    {
                        Value = 0.00001, //Does not show when set to 0
                        StrokeThickness = 1,
                        Stroke = new SolidColorBrush(System.Windows.Media.Color.FromRgb(64, 79, 86))
                    }

                }
            });

            dataTree.Dock = DockStyle.Fill;
            elutionPlot.Dock = DockStyle.Fill;
            elutionPlot.Hoverable = false;
            elutionPlot.DataTooltip = null;
            elutionPlot.DisableAnimations = true;
            spectrumPlot.Dock = DockStyle.Fill;
            spectrumPlot.Hoverable = false;
            spectrumPlot.DataTooltip = null;
            spectrumPlot.DisableAnimations = true;
            statsPlot1.Dock = DockStyle.Fill;
            statsPlot1.Hoverable = false;
            statsPlot1.DataTooltip = null;
            statsPlot2.Dock = DockStyle.Fill;
            statsPlot2.Hoverable = false;
            statsPlot2.DataTooltip = null;
            statsPlot3.Dock = DockStyle.Fill;
            statsPlot3.Hoverable = false;
            statsPlot3.DataTooltip = null;
            statsPlot3.DisableAnimations = true;
            statsPlot4.Dock = DockStyle.Fill;
            elutionPlot.AxisY.Clear();
            elutionPlot.AxisY.Add(new LiveCharts.Wpf.Axis
            {
                Title = "Intensity",
                MinValue = 0.0,
                FontSize = 15,
                Foreground = System.Windows.Media.Brushes.Black,
                Sections = new SectionsCollection
                {
                    new AxisSection
                    {
                        Value = 0.00001, //Does not show when set to 0
                        StrokeThickness = 1,
                        Stroke = new SolidColorBrush(System.Windows.Media.Color.FromRgb(64, 79, 86))
                    }

                }

            });
            elutionPlot.AxisX.Clear();
            elutionPlot.AxisX.Add(new LiveCharts.Wpf.Axis
            {
                Title = "Time (min.)",
                FontSize = 15,
                Foreground = System.Windows.Media.Brushes.Black,
                Sections = new SectionsCollection
                {
                    new AxisSection
                    {
                        Value = 0.00001, //Does not show when set to 0
                        StrokeThickness = 1,
                        Stroke = new SolidColorBrush(System.Windows.Media.Color.FromRgb(64, 79, 86))
                    }

                }
            });
            spectrumPlot.AxisY.Clear();
            spectrumPlot.AxisY.Add(new LiveCharts.Wpf.Axis
            {
                Title = "Intensity",
                MinValue = 0.0,
                FontSize = 15,
                Foreground = System.Windows.Media.Brushes.Black,
                Sections = new SectionsCollection
                {
                    new AxisSection
                    {
                        Value = 0.00001, //Does not show when set to 0
                        StrokeThickness = 1,
                        Stroke = new SolidColorBrush(System.Windows.Media.Color.FromRgb(64, 79, 86))
                    }

                }
            });
            spectrumPlot.AxisX.Clear();
            spectrumPlot.AxisX.Add(new LiveCharts.Wpf.Axis
            {
                Title = "m/z",
                MinValue = 100.0,
                FontSize = 15,
                Foreground = System.Windows.Media.Brushes.Black,
                Sections = new SectionsCollection
                {
                    new AxisSection
                    {
                        Value = 0.00001, //Does not show when set to 0
                        StrokeThickness = 1,
                        Stroke = new SolidColorBrush(System.Windows.Media.Color.FromRgb(64, 79, 86))
                    }

                }
                //MaxValue = 3000.0
            });

            statsPlot1.AxisX.Clear();
            statsPlot1.AxisX.Add(new LiveCharts.Wpf.Axis
            {
                Title = "Product Ion m/z",
                MinValue = 0,
                FontSize = 15,
                Foreground = System.Windows.Media.Brushes.Black,
                FontFamily = new System.Windows.Media.FontFamily("Calibri"),
                Sections = new SectionsCollection
                {
                    new AxisSection
                    {
                        Value = 0.00001, //Does not show when set to 0
                        StrokeThickness = 1,
                        Stroke = new SolidColorBrush(System.Windows.Media.Color.FromRgb(64, 79, 86))
                    }

                }
            });
            statsPlot1.AxisY.Clear();
            statsPlot1.AxisY.Add(new LiveCharts.Wpf.Axis
            {
                Title = "Count",
                MinValue = 0,
                FontSize = 15,
                Foreground = System.Windows.Media.Brushes.Black,
                Sections = new SectionsCollection
                {
                    new AxisSection
                    {
                        Value = 0.00001, //Does not show when set to 0
                        StrokeThickness = 1,
                        Stroke = new SolidColorBrush(System.Windows.Media.Color.FromRgb(64, 79, 86))
                    }

                }
            });

            statsPlot2.AxisX.Clear();
            statsPlot2.AxisX.Add(new LiveCharts.Wpf.Axis
            {
                Title = "Sequence Coverage (%)",
                MinValue = 0,
                MaxValue = 100,
                FontSize = 15,
                Foreground = System.Windows.Media.Brushes.Black,
                Sections = new SectionsCollection
                {
                    new AxisSection
                    {
                        Value = 0.00001, //Does not show when set to 0
                        StrokeThickness = 1,
                        Stroke = new SolidColorBrush(System.Windows.Media.Color.FromRgb(64, 79, 86))
                    }

                }
            });
            statsPlot2.AxisY.Clear();
            statsPlot2.AxisY.Add(new LiveCharts.Wpf.Axis
            {
                Title = "Count",
                MinValue = 0,
                FontSize = 15,
                Foreground = System.Windows.Media.Brushes.Black,
                Sections = new SectionsCollection
                {
                    new AxisSection
                    {
                        Value = 0.00001, //Does not show when set to 0
                        StrokeThickness = 1,
                        Stroke = new SolidColorBrush(System.Windows.Media.Color.FromRgb(64, 79, 86))
                    }

                }
            });

            statsPlot3.AxisX.Clear();
            statsPlot3.AxisX.Add(new LiveCharts.Wpf.Axis
            {
                Title = "Precursor m/z",
                FontSize = 15,
                Foreground = System.Windows.Media.Brushes.Black,
                Sections = new SectionsCollection
                {
                    new AxisSection
                    {
                        Value = 0.00001, //Does not show when set to 0
                        StrokeThickness = 1,
                        Stroke = new SolidColorBrush(System.Windows.Media.Color.FromRgb(64, 79, 86))
                    }

                }
            });
            statsPlot3.AxisY.Clear();
            statsPlot3.AxisY.Add(new LiveCharts.Wpf.Axis
            {
                Title = "Mass error (ppm)",
                MaxValue = 20,
                MinValue = -20,
                FontSize = 15,
                Foreground = System.Windows.Media.Brushes.Black,
                Sections = new SectionsCollection
                {
                    new AxisSection
                    {
                        Value = 0.00001, //Does not show when set to 0
                        StrokeThickness = 1,
                        Stroke = new SolidColorBrush(System.Windows.Media.Color.FromRgb(64, 79, 86))
                    }

                }
            });

            //resultViewBox.Items.AddRange(new string[] { "GlycoPSMs", "Glycopeptides", "Glycosites", "Glycoproteins", "Glycans" });
            dataGridView1.DataSource = bindingSource1;
            //resultViewBox.SelectedItem = resultViewBox.Items[0];

            //Testing data grid instead of data tree
            dataTree.Visible = false;
            viewPeptidesDataGrid.Dock = DockStyle.Fill;
            viewPeptidesDataGrid.DataSource = bindingSource2;
        }

        private void byrstLoad_Click(object sender, EventArgs e)
        {
            //Read in data and store in new SQLite DB'

            if (!dataFromDB)
            {
                BuildSQLiteResults();
            }
            else
            {
                fillTableNames();
                resultViewBox.SelectedItem = resultViewBox.Items[resultViewBox.Items.IndexOf("GlycoPSMs")];
                fillStatsComboBox();
                statsComboBox.SelectedItem = statsComboBox.Items[statsComboBox.Items.IndexOf("All Files")];
                fillResultTable("SELECT * FROM AllGlycoPSMs");
                fillViewPeptidesDataGrid("SELECT * FROM AllGlycoPSMs");          
                fillFragHistPlot("", true);
                fillSeqCoveragePlot("", true);
                fillPrecursorMassErrorPlot("", true);
                fillGlycanTypePieChart("", true);
                dataUploaded = true;
                fillInSourceFragData();
                AddProgressText("\nFinished reading data from database.");
            }
        }

        private void byrsltFiles_SelectedIndexChanged(object sender, EventArgs e)
        {

        }

        private void Form1_DragEnter(object sender, DragEventArgs e)
        {

        }

        private void Form1_DragDrop(object sender, DragEventArgs e)
        {

        }

        private void listView1_SelectedIndexChanged(object sender, EventArgs e)
        {

        }

        // Helper function for ByonicDataReader.Parse()
        // Reads byonic results from .byrslt files and stores data in SQLite Table
        private void BuildSQLiteResults()
        {
            // files: { FilePath, Condition, Replicate, isControl }
            var files = new List<string[]>();

            for (int i = 0; i < dataUpload.Rows.Count; i++)
            {
                // Add checks to make sure form is filled out
                string[] file = { dataUpload.Rows[i].Cells[0].Value.ToString(),
                                  dataUpload.Rows[i].Cells[1].Value.ToString(),
                                  dataUpload.Rows[i].Cells[2].Value.ToString(),
                                  dataUpload.Rows[i].Cells[3].Value.ToString(),
                                  (dataUpload.Rows[i].Cells[4].Value != System.DBNull.Value).ToString()
                                };

                files.Add(file);
            }

            // Perform SQLite Table Build in Separate Thread
            AddProgressText("Starting data upload...");
            var bynReader = new ByonicDataReader(files, outputPath.Text);
            bynReader.scoreFilter = (double)scoreFilter.Value;
            bynReader.deltaModFilter = (double)deltaModScoreFilter.Value;
            bynReader.logProbFilter = (double)logProbFilter.Value;
            bynReader.pepLengthFilter = (int)pepLengthFilter.Value;
            bynReader.glycanCountFilter = Int32.Parse(modCountFilter.Text);
            bynReader.organism = organismText.Text;
            bynReader.performProteinInference = PerformProteinInference.Checked;
            bynReader.performGlycanLocalization = PerformGlycanLocalization.Checked;
            bynReader.IdentifyForInsourceFragments = IdentifyInsourceFragments.Checked;

            if (PerformProteinInference.Checked)
            {
                bynReader.AdditionalPeptidesPath = AdditionalPeptidesPath.Text;
                bynReader.fastaFile = fastaPath.Text;
                bynReader.proteaseString = ProteaseCB.Items[ProteaseCB.SelectedIndex].ToString();
                bynReader.MaxMissedCleavage = Int32.Parse(miscleaveCB.Items[miscleaveCB.SelectedIndex].ToString());
            }

            bynReader.UpdateProgress += HandleUpdateProgress;
            bynReader.finish += HandleFinishedUpload;
            Thread thread = new Thread(bynReader.getNewData);
            thread.Start();
            ChangeButtonStatus(false);
        }

        private void AssociateRawFile(string path)
        {
            var file = Path.GetFileNameWithoutExtension(path);

            for (int i = 0; i < dataUpload.Rows.Count; i++)
            {
                if (dataUpload.Rows[i].Cells[0].Value.ToString().Contains(file + ".raw") || dataUpload.Rows[i].Cells[0].Value.ToString().Contains(file + ".mgf"))
                {
                    dataUpload.Rows[i].Cells[1].Value = path;
                }
            }
        }

        private void GlycoTools_Load(object sender, EventArgs e)
        {
            //CreateGraph(plotBox);
            //SetSize();
        }

        private void CreateGraph(ZedGraphControl zgc)
        {

        }

        private void GlycoTools_Resize(object sender, EventArgs e)
        {
            //SetSize();
        }

        private void SetSize()
        {

            //plotBox.Location = new Point(10, 10);
            //plotBox.Size = new Size(ClientRectangle.Width - 20, ClientRectangle.Height - 20);
        }

        private void HandleUpdateProgress(object sender, ProgressEventArgs e)
        {
            AddProgressText(e.progressText);
        }

        private void HandleFinishedUpload(object sender, DataReturnArgs e)
        {

            //Consider putting these into a task
            AddExperimenDesignTableToSQLDB();
            fillTableNames();
            resultViewBox.SelectedItem = resultViewBox.Items[resultViewBox.Items.IndexOf("GlycoPSMs")];
            fillStatsComboBox();
            statsComboBox.SelectedItem = statsComboBox.Items[statsComboBox.Items.IndexOf("All Files")];
            fillResultTable("SELECT * FROM AllGlycoPSMs");
            fillViewPeptidesDataGrid("SELECT * FROM AllGlycoPSMs");
            BatchAnnotateSpectra();
            AddProgressText("\nFinished Reading Data.");
            fillFragHistPlot("", true);
            fillSeqCoveragePlot("", true);
            fillPrecursorMassErrorPlot("", true);
            fillGlycanTypePieChart("", true);
            dataUploaded = true;
            if (IdentifyInsourceFragments.Checked)
            {
                fillInSourceFragData();
            }           

            ChangeButtonStatus(true);
        }

        private void AddExperimenDesignTableToSQLDB()
        {
            var connection = new SQLiteConnection(string.Format("Data Source={0}; Version=3;", outputPath.Text + "\\MyDatabase.sqlite"));
            connection.Open();
            using (var transaction1 = connection.BeginTransaction())
            {
                var commandString = "CREATE TABLE IF NOT EXISTS ExperimentDesign (Byrslt STRING, Raw STRING, Condition STRING, Replicate STRING, IsControl BOOL)";
                var command = new SQLiteCommand(commandString, connection);
                command.ExecuteNonQuery();

                foreach (DataGridViewRow row in dataUpload.Rows)
                {
                    var byrslt = row.Cells[0].Value.ToString();
                    var raw = row.Cells[1].Value.ToString();
                    var condition = row.Cells[2].Value.ToString();
                    var replicate = row.Cells[3].Value.ToString();
                    var isControl = Convert.ToBoolean(row.Cells[4].Value);

                    var insertString = string.Format("INSERT INTO ExperimentDesign('Byrslt', 'Raw', 'Condition', 'Replicate', 'IsControl') VALUES ('{0}', '{1}', '{2}', '{3}', '{4}')", byrslt, raw, condition, replicate, isControl);
                    var insertCommand = new SQLiteCommand(insertString, connection);
                    var reader = insertCommand.ExecuteReader();
                }
                transaction1.Commit();
            }
            
            connection.Close();

        }

        private void fillInSourceFragData()
        {
            var connection = new SQLiteConnection(string.Format("Data Source={0}; Version=3;", outputPath.Text + "\\MyDatabase.sqlite"));

            var selectCommand = "SELECT * FROM InSourceFragData";
            dataAdapter = new SQLiteDataAdapter(selectCommand, connection);

            DataTable table = new DataTable
            {
                Locale = CultureInfo.InvariantCulture
            };

            dataAdapter.Fill(table);
            inSourceFragBindingSource.DataSource = table;
            inSourceFragData.AutoResizeColumns(DataGridViewAutoSizeColumnsMode.AllCells);
            connection.Close();
        }

        private void fillFragHistPlot(string file, bool allFiles)
        {
            if (InvokeRequired)
            {
                statsPlot1.Invoke(new Action<string, bool>(fillFragHistPlot), file, allFiles);
                return;
            }

            while (statsPlot1.Series.Count > 0)
                statsPlot1.Series.RemoveAt(statsPlot1.Series.Count() - 1);

            var connection = new SQLiteConnection("Data Source=" + outputPath.Text + "\\MyDatabase.sqlite; Version=3;");
            connection.Open();

            //Get min
            var min = 0.0;
            var minCommandString = "";
            if (allFiles)
            {
                minCommandString = "SELECT * FROM FragmentIons ORDER BY MZ ASC";
            }
            else
            {
                minCommandString = string.Format("SELECT * FROM FragmentIons WHERE FragmentIons.File='{0}' ORDER BY MZ ASC", file);
            }

            var minCommand = new SQLiteCommand(minCommandString, connection);
            var minReader = minCommand.ExecuteReader();

            while (minReader.Read())
            {
                min = double.Parse(minReader["MZ"].ToString());
                break;
            }
            minReader.Close();

            //Get max
            var max = 0.0;
            var maxCommandString = "";
            if (allFiles)
            {
                maxCommandString = "SELECT * FROM FragmentIons ORDER BY MZ DESC";
            }
            else
            {
                maxCommandString = string.Format("SELECT * FROM FragmentIons WHERE FragmentIons.File='{0}' ORDER BY MZ DESC", file);
            }

            var maxCommand = new SQLiteCommand(maxCommandString, connection);
            var maxReader = maxCommand.ExecuteReader();

            while (maxReader.Read())
            {
                max = double.Parse(maxReader["MZ"].ToString());
                break;
            }
            maxReader.Close();

            //pepHistReader
            var pephistCommandString = "";
            if (allFiles)
            {
                pephistCommandString = "SELECT * FROM FragmentIons WHERE FragmentIons.IonType='Backbone'";
            }
            else
            {
                pephistCommandString = string.Format("SELECT * FROM FragmentIons WHERE FragmentIons.File='{0}' AND FragmentIons.IonType='Backbone'", file);
            }
            var pepHistCommand = new SQLiteCommand(pephistCommandString, connection);
            var pepHistReader = pepHistCommand.ExecuteReader();
            var pepHist = GetBinCountsFromSqlReader(30, pepHistReader, min, max, "MZ");
            pepHistReader.Close();

            //yIonHistReader
            var yIonHistCommandString = "";
            if (allFiles)
            {
                yIonHistCommandString = "SELECT * FROM FragmentIons WHERE FragmentIons.IonType='YIon'";
            }
            else
            {
                yIonHistCommandString = string.Format("SELECT * FROM FragmentIons WHERE FragmentIons.File='{0}' AND FragmentIons.IonType='YIon'", file);
            }
            var yIonHistCommand = new SQLiteCommand(yIonHistCommandString, connection);
            var yIonHistReader = yIonHistCommand.ExecuteReader();
            var yIonHist = GetBinCountsFromSqlReader(30, yIonHistReader, min, max, "MZ");
            yIonHistReader.Close();

            //OxoniumHistReader
            var oxonHistCommandString = "";
            if (allFiles)
            {
                oxonHistCommandString = "SELECT * FROM FragmentIons WHERE FragmentIons.IonType='Oxonium'";
            }
            else
            {
                oxonHistCommandString = string.Format("SELECT * FROM FragmentIons WHERE FragmentIons.File='{0}' AND FragmentIons.IonType='Oxonium'", file);
            }
            var oxonHistCommand = new SQLiteCommand(oxonHistCommandString, connection);
            var oxonHistReader = oxonHistCommand.ExecuteReader();
            var oxonHist = GetBinCountsFromSqlReader(30, oxonHistReader, min, max, "MZ");
            oxonHistReader.Close();
            connection.Close();
            GC.Collect();

            var pepData = new ChartValues<ObservablePoint>();
            for (int i = 0; i < pepHist[0].Count(); i++)
            {
                pepData.Add(new ObservablePoint
                {
                    X = pepHist[0][i],
                    Y = pepHist[1][i] // pepHist[1].Sum()
                });
            }

            statsPlot1.Series.Add(new ColumnSeries
            {
                Title = "Backbone Fragment Ions",
                Values = pepData,
                ColumnPadding = -7,
                Fill = System.Windows.Media.Brushes.LightCoral,
                Stroke = System.Windows.Media.Brushes.White,
                StrokeThickness = 0.5,
                //Width = 1000,

            });

            var yIonData = new ChartValues<ObservablePoint>();
            for (int i = 0; i < yIonHist[0].Count(); i++)
            {
                yIonData.Add(new ObservablePoint
                {
                    X = yIonHist[0][i],
                    Y = yIonHist[1][i] // yIonHist[1].Sum()
                });
            }

            statsPlot1.Series.Add(new ColumnSeries
            {
                Title = "Y Ions",
                Values = yIonData,
                ColumnPadding = -7,
                Fill = System.Windows.Media.Brushes.DodgerBlue,
                Stroke = System.Windows.Media.Brushes.White,
                StrokeThickness = 0.5,
                //Width = 1000,

            });

            var oxoniumIonData = new ChartValues<ObservablePoint>();
            for (int i = 0; i < oxonHist[0].Count(); i++)
            {
                oxoniumIonData.Add(new ObservablePoint
                {
                    X = oxonHist[0][i],
                    Y = oxonHist[1][i] // oxonHist[1].Sum()
                });
            }

            statsPlot1.Series.Add(new ColumnSeries
            {
                Title = "Oxonium Ions",
                Values = oxoniumIonData,
                ColumnPadding = -7,
                Fill = System.Windows.Media.Brushes.SeaGreen,
                Stroke = System.Windows.Media.Brushes.White,
                StrokeThickness = 0.5,
                //Width = 1000,               

            });

            statsPlot1.DefaultLegend.FontSize = 15;
            statsPlot1.DefaultLegend.FontFamily = new System.Windows.Media.FontFamily("Calibri");
            statsPlot1.DefaultLegend.BorderBrush = System.Windows.Media.Brushes.Gray;
            statsPlot1.LegendLocation = LiveCharts.LegendLocation.Bottom;


        }

        private void fillSeqCoveragePlot(string file, bool allFiles)
        {
            if (InvokeRequired)
            {
                statsPlot2.Invoke(new Action<string, bool>(fillSeqCoveragePlot), file, allFiles);
                return;
            }

            while (statsPlot2.Series.Count > 0)
                statsPlot2.Series.RemoveAt(statsPlot2.Series.Count() - 1);

            var connection = new SQLiteConnection("Data Source=" + outputPath.Text + "\\MyDatabase.sqlite; Version=3;");
            connection.Open();

            //Pep Reader
            var pepCommandString = "";
            if (allFiles)
            {
                pepCommandString = "SELECT * FROM AnnotatedPeptides";
            }
            else
            {
                pepCommandString = string.Format("SELECT * FROM AnnotatedPeptides WHERE AnnotatedPeptides.File='{0}'", file);
            }
            var pepCommand = new SQLiteCommand(pepCommandString, connection);
            var pepReader = pepCommand.ExecuteReader();
            var pepHistData = GetBinCountsFromSqlReader(30, pepReader, 0, 1, "SeqCoverage");
            pepReader.Close();

            //Glycan Reader
            var glyCommandString = "";
            if (allFiles)
            {
                glyCommandString = "SELECT * FROM AnnotatedPeptides";
            }
            else
            {
                glyCommandString = string.Format("SELECT * FROM AnnotatedPeptides WHERE AnnotatedPeptides.File='{0}'", file);
            }
            var glyCommand = new SQLiteCommand(glyCommandString, connection);
            var glyReader = glyCommand.ExecuteReader();
            var glycanHistData = GetBinCountsFromSqlReader(30, glyReader, 0, 1, "GlycanSeqCoverage");
            glyReader.Close();

            var pepData = new ChartValues<ObservablePoint>();
            for (int i = 0; i < pepHistData[0].Count(); i++)
            {
                pepData.Add(new ObservablePoint
                {
                    X = pepHistData[0][i] * 100,
                    Y = pepHistData[1][i]// pepHist[1].Sum()
                });
            }

            statsPlot2.Series.Add(new ColumnSeries
            {
                Title = "Peptide Backbone",
                Values = pepData,
                ColumnPadding = -5,
                Fill = System.Windows.Media.Brushes.LightCoral,
                Stroke = System.Windows.Media.Brushes.White,
                StrokeThickness = 1,
                //Width = 1000,
            });

            var glyData = new ChartValues<ObservablePoint>();
            for (int i = 0; i < glycanHistData[0].Count(); i++)
            {
                glyData.Add(new ObservablePoint
                {
                    X = glycanHistData[0][i] * 100,
                    Y = glycanHistData[1][i] // pepHist[1].Sum()
                });
            }

            statsPlot2.Series.Add(new ColumnSeries
            {
                Title = "Glycan",
                Values = glyData,
                ColumnPadding = -5,
                Fill = System.Windows.Media.Brushes.DodgerBlue,
                Stroke = System.Windows.Media.Brushes.White,
                StrokeThickness = 1,
                //Width = 1000,
            });

            statsPlot2.DefaultLegend.FontSize = 15;
            statsPlot2.LegendLocation = LiveCharts.LegendLocation.Bottom;

            connection.Close();
            GC.Collect();
        }

        private void fillGlycanTypePieChart(string file, bool allFiles)
        {
            if (InvokeRequired)
            {
                statsPlot4.Invoke(new Action<string, bool>(fillGlycanTypePieChart), file, allFiles);
                return;
            }

            while (statsPlot4.Series.Count > 0)
                statsPlot4.Series.RemoveAt(statsPlot4.Series.Count() - 1);

            var connection = new SQLiteConnection("Data Source=" + outputPath.Text + "\\MyDatabase.sqlite; Version=3;");
            connection.Open();

            var commandString = "";
            if (allFiles)
            {
                commandString = "SELECT * FROM AllGlycoPSMs";
            }
            else
            {
                commandString = string.Format("SELECT * FROM AllGlycoPSMs WHERE AllGlycoPSMs.File='{0}'", file);
            }

            var command = new SQLiteCommand(commandString, connection);
            var reader = command.ExecuteReader();

            var data = new Dictionary<string, int>();
            while (reader.Read())
            {
                var types = reader["GlycanTypes"].ToString().Split(';');

                foreach (var type in types)
                {
                    if (data.ContainsKey(type))
                    {
                        data[type]++;
                    }
                    else
                    {
                        data.Add(type, 1);
                    }
                }
            }

            statsPlot4.Series = new SeriesCollection { };

            Func<ChartPoint, string> labelPoint = chartPoint => string.Format("{0} ({1:P})", chartPoint.Y, chartPoint.Participation);

            foreach (var pair in data)
            {
                switch (pair.Key)
                {
                    case "High Mannose":
                        statsPlot4.Series.Add(new PieSeries
                        {
                            Title = pair.Key,
                            Values = new ChartValues<double> { pair.Value },
                            DataLabels = true,
                            LabelPoint = labelPoint,
                            Fill = new SolidColorBrush
                            {
                                Color = System.Windows.Media.Color.FromRgb(46, 139, 87)
                            },
                        });
                        break;

                    case "Sialylated":
                        statsPlot4.Series.Add(new PieSeries
                        {
                            Title = pair.Key,
                            Values = new ChartValues<double> { pair.Value },
                            DataLabels = true,
                            LabelPoint = labelPoint,
                            Fill = new SolidColorBrush
                            {
                                Color = System.Windows.Media.Color.FromRgb(0, 90, 156)
                            },
                        });
                        break;

                    case "Fucosylated":
                        statsPlot4.Series.Add(new PieSeries
                        {
                            Title = pair.Key,
                            Values = new ChartValues<double> { pair.Value },
                            DataLabels = true,
                            LabelPoint = labelPoint,
                            Fill = new SolidColorBrush
                            {
                                Color = System.Windows.Media.Color.FromRgb(240, 128, 128)
                            },
                        });
                        break;

                    case "Complex/Hybrid":
                        statsPlot4.Series.Add(new PieSeries
                        {
                            Title = pair.Key,
                            Values = new ChartValues<double> { pair.Value },
                            DataLabels = true,
                            LabelPoint = labelPoint,
                            Fill = new SolidColorBrush
                            {
                                Color = System.Windows.Media.Color.FromRgb(211, 211, 211)
                            },
                        });
                        break;

                    case "Paucimannose":
                        statsPlot4.Series.Add(new PieSeries
                        {
                            Title = pair.Key,
                            Values = new ChartValues<double> { pair.Value },
                            DataLabels = true,
                            LabelPoint = labelPoint,
                            Fill = new SolidColorBrush
                            {
                                Color = System.Windows.Media.Color.FromRgb(152, 251, 152)
                            },
                        });
                        break;

                    case "Phosphomannose":
                        statsPlot4.Series.Add(new PieSeries
                        {
                            Title = pair.Key,
                            Values = new ChartValues<double> { pair.Value },
                            DataLabels = true,
                            LabelPoint = labelPoint,
                            Fill = new SolidColorBrush
                            {
                                Color = System.Windows.Media.Color.FromRgb(255, 99, 71)
                            },
                        });
                        break;
                }



            }

            statsPlot4.DefaultLegend.FontSize = 15;
            statsPlot4.LegendLocation = LiveCharts.LegendLocation.Bottom;
        }

        private void fillPrecursorMassErrorPlot(string file, bool allFiles)
        {
            if (InvokeRequired)
            {
                statsPlot3.Invoke(new Action<string, bool>(fillPrecursorMassErrorPlot), file, allFiles);
                return;
            }

            while (statsPlot3.Series.Count > 0)
                statsPlot3.Series.RemoveAt(statsPlot3.Series.Count() - 1);

            var connection = new SQLiteConnection("Data Source=" + outputPath.Text + "\\MyDatabase.sqlite; Version=3;");
            connection.Open();

            var commandString = "";
            if (allFiles)
            {
                commandString = "SELECT * FROM AllGlycoPSMs";
            }
            else
            {
                commandString = string.Format("SELECT * FROM AllGlycoPSMs WHERE AllGlycoPSMs.File='{0}'", file);
            }

            var command = new SQLiteCommand(commandString, connection);
            var reader = command.ExecuteReader();

            var pepData = new ChartValues<ObservablePoint>();
            while (reader.Read())
            {
                var error = double.Parse(reader["ppmerr"].ToString());
                var mz = double.Parse(reader["ObsMZ"].ToString());

                pepData.Add(new ObservablePoint
                {
                    X = mz,
                    Y = error
                });
            }

            statsPlot3.Series.Add(new ScatterSeries
            {
                Values = pepData,
                MaxPointShapeDiameter = 10,
                Fill = new SolidColorBrush
                {
                    Color = System.Windows.Media.Color.FromRgb(0, 90, 156),
                    Opacity = .4
                },
                Stroke = System.Windows.Media.Brushes.DodgerBlue,
                StrokeThickness = 1
            });

            connection.Close();
            GC.Collect();

        }

        private void fillTableNames()
        {
            resultViewBox.Items.AddRange(new string[] { "GlycoPSMs", "Glycopeptides", "Glycoproteins", "Glycosites", "Glycans" });
            resultViewBox.SelectedIndex = 0;
            conditionComboBox.Items.Add("All");
            conditionComboBox.SelectedIndex = 0;
            //replicateComboBox.Items.Add("All");
            //replicateComboBox.SelectedIndex = 0;

            var connection = new SQLiteConnection("Data Source=" + outputPath.Text + "\\MyDatabase.sqlite; Version=3;");
            connection.Open();

            var command = new SQLiteCommand("SELECT DISTINCT Condition FROM AllGlycoPSMs", connection);
            var reader = command.ExecuteReader();

            while (reader.Read())
            {
                var condition = reader["Condition"].ToString();
                conditionComboBox.Items.Add(condition);
            }

            connection.Close();
            GC.Collect();
        }

        private void addResultTable(string[] names)
        {
            if (InvokeRequired)
            {
                resultViewBox.Invoke(new Action<string[]>(addResultTable), names);
                return;
            }


        }

        private void fillStatsComboBox()
        {
            statsComboBox.Items.Add("All Files");

            var connection = new SQLiteConnection("Data Source=" + outputPath.Text + "\\MyDatabase.sqlite; Version=3;");
            connection.Open();

            var command = new SQLiteCommand("SELECT DISTINCT File FROM AllGlycoPSMS", connection);
            var reader = command.ExecuteReader();

            while (reader.Read())
            {
                var file = reader["File"];
                statsComboBox.Items.Add(file);
            }

            connection.Close();
        }

        private void fillViewPeptidesDataGrid(string selectCommand)
        {
            if (InvokeRequired)
            {
                dataGridView1.Invoke(new Action<string>(fillViewPeptidesDataGrid), selectCommand);
                return;
            }

            var connection = new SQLiteConnection("Data Source=" + outputPath.Text + "\\MyDatabase.sqlite; Version=3;");
            dataAdapter = new SQLiteDataAdapter(selectCommand, connection);

            DataTable table = new DataTable
            {
                Locale = CultureInfo.InvariantCulture
            };

            dataAdapter.Fill(table);
            bindingSource2.DataSource = table;
            dataGridView1.AutoResizeColumns(DataGridViewAutoSizeColumnsMode.AllCells);

            foreach (DataGridViewColumn column in viewPeptidesDataGrid.Columns)
            {
                column.Width = 100;
            }

            connection.Close();
            GC.Collect();

        }

        private void fillResultTable(string selectCommand)
        {
            var connection = new SQLiteConnection("Data Source=" + outputPath.Text + "\\MyDatabase.sqlite; Version=3;");
            dataAdapter = new SQLiteDataAdapter(selectCommand, connection);

            //SQLiteCommandBuilder commandBuilder = new SQLiteCommandBuilder(dataAdapter);

            DataTable table = new DataTable
            {
                Locale = CultureInfo.InvariantCulture
            };

            dataAdapter.Fill(table);
            bindingSource1.DataSource = table;
            dataGridView1.AutoResizeColumns(DataGridViewAutoSizeColumnsMode.AllCells);

            foreach (DataGridViewColumn column in dataGridView1.Columns)
            {
                column.Width = 100;
            }

            connection.Close();
            GC.Collect();

        }

        private void ChangeButtonStatus(bool enabled)
        {
            if (InvokeRequired)
            {
                byrsltLoad.Invoke(new Action<bool>(ChangeButtonStatus), enabled);
                return;
            }
            byrsltLoad.Enabled = enabled;
        }

        private void AddProgressText(string progressValue)
        {
            if (InvokeRequired)
            {
                progressText.Invoke(new Action<string>(AddProgressText), progressValue);
                return;
            }
            progressText.AppendText(string.Format("\n{0}", progressValue));
            progressText.AppendText(Environment.NewLine);
        }

        private void dataUpload_CellContentClick(object sender, DataGridViewCellEventArgs e)
        {

        }

        private void dataUpload_DragLeave(object sender, EventArgs e)
        {

        }

        private void dataUpload_DragDrop(object sender, DragEventArgs e)
        {
            if (!e.Data.GetDataPresent(DataFormats.FileDrop))
                return;

            string[] files = e.Data.GetData(DataFormats.FileDrop) as string[];
            if (files == null)
                return;

            foreach (string file in files.Where(f => Path.GetExtension(f).Equals(".byrslt")))
            {
                string[] row = { file };

                dataUpload.Rows.Add(row);

                if (String.IsNullOrEmpty(outputPath.Text))
                {
                    outputPath.Text = Path.GetDirectoryName(file);
                }

                textBox1.Visible = false;
                textBox2.Visible = false;
                textBox3.Visible = false;
                //textBox4.Visible = false;
                //textBox5.Visible = false;
            }

            foreach (string file in files.Where(f => Path.GetExtension(f).Equals(".raw")))
            {
                AssociateRawFile(file);
            }

            foreach(string file in files.Where(f => Path.GetExtension(f).Equals(".sqlite")))
            {
                try
                {
                    ExperimentDesignFromSQLDB(file);

                    textBox1.Visible = false;
                    textBox2.Visible = false;
                    textBox3.Visible = false;
                    //textBox4.Visible = false;
                    //textBox5.Visible = false;

                    dataFromDB = true;
                    dataUpload.ReadOnly = true;

                    if (String.IsNullOrEmpty(outputPath.Text))
                    {
                        outputPath.Text = Path.GetDirectoryName(file);
                    }
                }
                catch(Exception exc)
                {
                    AddProgressText("Failed to read data from " + Path.GetFileName(file));
                }                
            }
        }

        private void ExperimentDesignFromSQLDB(string file)
        {
            if (InvokeRequired)
            {
                dataUpload.Invoke(new Action<string>(ExperimentDesignFromSQLDB), file);
                return;
            }

            try
            {
                dataUpload.DataSource = dataUploadBindingSouce;
                dataUpload.Columns.Clear();
                dataUpload.Refresh();

                var connection = new SQLiteConnection(string.Format("Data Source={0}; Version=3;", file));

                var selectCommand = "SELECT Byrslt as 'Byonic Results (.byrslt)', Raw as 'Raw Files', Condition as 'Condition', Replicate as 'Replicate', IsControl as 'Control' FROM ExperimentDesign";

                dataAdapter = new SQLiteDataAdapter(selectCommand, connection);

                DataTable table = new DataTable
                {
                    Locale = CultureInfo.InvariantCulture
                };

                dataAdapter.Fill(table);
                dataUploadBindingSouce.DataSource = null;
                dataUploadBindingSouce.DataSource = table;
                connection.Close();
            }
            catch (Exception e)
            {
                //AddProgressText("Failed to get data from database.");
            }
        }

        private void dataUpload_DragEnter(object sender, DragEventArgs e)
        {
            e.Effect = DragDropEffects.None;

            if (!e.Data.GetDataPresent(DataFormats.FileDrop))
                return;

            string[] files = e.Data.GetData(DataFormats.FileDrop) as string[];
            if (files == null)
                return;

            if (files.Any(f => Path.GetExtension(f).Equals(".byrslt")))
                e.Effect = DragDropEffects.Link;

            if (files.Any(f => Path.GetExtension(f).Equals(".raw")))
                e.Effect = DragDropEffects.Link;

            if (files.Any(f => Path.GetExtension(f).Equals(".sqlite")))
                e.Effect = DragDropEffects.Link;
        }

        private void populateTreeView(Dictionary<string, Dictionary<string, List<string>>> data)
        {
            if (InvokeRequired)
            {
                dataTree.Invoke(new Action<Dictionary<string, Dictionary<string, List<string>>>>(populateTreeView), data);
                return;
            }
            dataTree.Nodes.Clear();

            int i = 0;

            //For each file
            foreach (KeyValuePair<string, Dictionary<string, List<string>>> file in data)
            {
                dataTree.Nodes.Add(file.Key);

                int j = 0;

                // For each protein
                foreach (KeyValuePair<string, List<string>> pair in file.Value)
                {
                    dataTree.Nodes[i].Nodes.Add(pair.Key);

                    //For each peptides
                    for (int k = 0; k < pair.Value.Count(); k++)
                    {
                        dataTree.Nodes[i].Nodes[j].Nodes.Add(pair.Value[k]);
                    }
                    j++;
                }
                i++;
            }
        }

        private void dataTree_NodeMouseClick(object sender, TreeNodeMouseClickEventArgs e)
        {
            if (e.Node.Level == 2)
            {
                var sequence = e.Node.Text;
                var prot = e.Node.Parent.Text;
                var file = e.Node.Parent.Parent.Text;
                //getPlotData(sequence, prot, file);
            }
        }

        private void getPlotDataFromSQL(int scanNum, string file)
        {
            if (InvokeRequired)
            {
                spectrumPlot.Invoke(new Action<int, string>(getPlotDataFromSQL), scanNum, file);
                return;
            }
            while (elutionPlot.Series.Count > 0)
                elutionPlot.Series.RemoveAt(elutionPlot.Series.Count() - 1);

            while (spectrumPlot.Series.Count > 0)
                spectrumPlot.Series.RemoveAt(spectrumPlot.Series.Count() - 1);

            while (spectrumPlot.VisualElements.Count() > 0)
                spectrumPlot.VisualElements.RemoveAt(spectrumPlot.VisualElements.Count() - 1);

            var sqlPath = string.Format("Data Source={0}\\MyDatabase.sqlite; Version=3;", outputPath.Text);
            SQLiteConnection sqlReader = new SQLiteConnection(sqlPath);
            sqlReader.Open();

            var query = string.Format("SELECT * FROM AllGlycoPSMs WHERE AllGlycoPSMs.ScanNum='{0}' AND AllGlycoPSMs.File='{1}' ORDER BY AllGlycoPSMs.Score DESC", scanNum, file);
            var command = new SQLiteCommand(query, sqlReader);
            var reader = command.ExecuteReader();

            var pepSequence = "";
            var glycans = new List<string>();
            var mods = new List<string>();
            var charge = 0;
            var RTs = new List<string>();
            var elutionIntensities = new List<string>();
            var glycoPSM = new GlycoPSM();

            while (reader.Read())
            {
                glycoPSM = SQLiteStringToGlycoPSM(reader);
                pepSequence = reader["Peptide"].ToString();
                glycans = reader["Glycans"].ToString().Split(';').ToList();
                mods = reader["ModsVar"].ToString().Split(';').ToList();
                charge = Int32.Parse(reader["Charge"].ToString());
                RTs = reader["RetentionTimes"].ToString().Split(',').ToList();
                elutionIntensities = reader["Intensities"].ToString().Split(',').ToList();
                break;
            }
            sqlReader.Close();
            GC.Collect();
            elutionPlot.Series = ChangeElutionPlot(RTs, elutionIntensities, file, scanNum);
            updateSpectrumPlotFromSQL(scanNum, file);

        }

        private void updateSpectrumPlotFromSQL(int scanNum, string file)
        {
            var rawFilePath = "";

            //var fileNoExt = Path.GetFileNameWithoutExtension(file);

            foreach (DataGridViewRow row in dataUpload.Rows)
            {
                var rowNoExt = Path.GetFileNameWithoutExtension(row.Cells[0].Value.ToString());
                if (rowNoExt.Equals(file))
                {
                    rawFilePath = row.Cells[1].Value.ToString();
                }
            }

            if (rawFilePath.Equals(""))
                return;

            var raw = new ThermoRawFile(rawFilePath);
            raw.Open();
            var spectrum = raw.GetSpectrum(scanNum);
            raw.Dispose();

            var sqlPath = string.Format("Data Source={0}\\MyDatabase.sqlite; Version=3;", outputPath.Text);
            SQLiteConnection sqlReader = new SQLiteConnection(sqlPath);
            sqlReader.Open();

            var query = string.Format("SELECT * FROM FragmentIons WHERE FragmentIons.ScanNumber='{0}' AND FragmentIons.File='{1}'", scanNum, file);
            var command = new SQLiteCommand(query, sqlReader);
            var reader = command.ExecuteReader();

            var pepFrags = new List<FragmentMatch>();
            var yIons = new List<FragmentMatch>();
            var oxoniumIons = new List<FragmentMatch>();

            while (reader.Read())
            {
                var fragmentName = reader["Fragment"].ToString();
                var fragmentType = reader["IonType"].ToString();
                var fragmentNumber = Int32.Parse(reader["FragmentNumber"].ToString());
                var fragmentMZ = double.Parse(reader["MZ"].ToString());
                var fragmentSignal = double.Parse(reader["Intensity"].ToString());
                var fragmentCharge = Int32.Parse(reader["Charge"].ToString());

                var fragMatch = new FragmentMatch(fragmentName, fragmentType, fragmentNumber, fragmentCharge, fragmentMZ, fragmentSignal);

                switch (fragmentType)
                {
                    case "Backbone":
                        pepFrags.Add(fragMatch);
                        break;

                    case "YIon":
                        yIons.Add(fragMatch);
                        break;

                    case "Oxonium":
                        oxoniumIons.Add(fragMatch);
                        break;
                }
            }

            sqlReader.Close();
            GC.Collect();

            var spectrumData = new List<List<FragmentMatch>>() { pepFrags, yIons, oxoniumIons };
            drawSpectrum(spectrum, spectrumData);


        }

        private void drawSpectrum(ThermoSpectrum spectrum, List<List<FragmentMatch>> spectrumData)
        {
            if (InvokeRequired)
            {
                spectrumPlot.Invoke(new Action<ThermoSpectrum, List<List<FragmentMatch>>>(drawSpectrum), spectrum, spectrumData);
                return;
            }

            // Full Spectrum
            var mzs = spectrum.GetMasses();
            var intensities = spectrum.GetIntensities();
            spectrumPlot.Series = new SeriesCollection();

            var plotData = new ChartValues<ObservablePoint>();

            for (int i = 0; i < mzs.Count(); i++)
            {
                plotData.Add(new ObservablePoint
                {
                    X = mzs[i],
                    Y = intensities[i],
                });
            }

            spectrumPlot.Series.Add(new ColumnSeries
            {
                Values = plotData,
                ColumnPadding = -0.5,
                Stroke = System.Windows.Media.Brushes.Black,
                StrokeThickness = 1000,
                Width = 100,

            });

            // Annotated Fragments
            var brushes = new List<System.Windows.Media.Brush>() { System.Windows.Media.Brushes.Red,
                                                                   System.Windows.Media.Brushes.Blue,
                                                                   System.Windows.Media.Brushes.Green,
                                                                   System.Windows.Media.Brushes.Brown,
                                                                   System.Windows.Media.Brushes.DarkGray};
            var j = 0;
            foreach (var spec in spectrumData)
            {

                var specIntensities = spec.Select(x => x.fragmentSignal).ToList();
                var specMZs = spec.Select(x => x.fragmentMZ).ToList();
                var labels = spec.Select(x => x.label).ToList();

                var specData = new ChartValues<ObservablePoint>();
                for (int i = 0; i < specMZs.Count(); i++)
                {
                    specData.Add(new ObservablePoint
                    {
                        X = specMZs[i],
                        Y = specIntensities[i]
                    });

                    var label = new VisualElement
                    {
                        X = specMZs[i] + 10,
                        Y = specIntensities[i],
                        VerticalAlignment = System.Windows.VerticalAlignment.Top,
                        HorizontalAlignment = System.Windows.HorizontalAlignment.Right,
                        UIElement = new TextBlock { Text = labels[i], RenderTransformOrigin = new System.Windows.Point(0.0, 0.9), RenderTransform = new RotateTransform(-90) }, // , LayoutTransform = new RotateTransform(-90)},                        
                        //UIElement.RenderTransform = new RotateTransform(-90),
                        //LayoutTransform = new RotateTransform(-90),
                    };
                    spectrumPlot.VisualElements.Add(label);
                }

                var title = "";
                switch (j)
                {
                    case 0:
                        title = "Peptide Backbone";
                        break;
                    case 1:
                        title = "Y Ions";
                        break;
                    case 2:
                        title = "Oxonium Ions";
                        break;
                }
                spectrumPlot.Series.Add(new ColumnSeries
                {
                    Title = title,
                    Values = specData,
                    ColumnPadding = -2,
                    Stroke = brushes[j],
                    StrokeThickness = 1000,
                    Width = 1000,

                });

                j++;
            }

            spectrumPlot.AxisY[0].MaxValue = spectrumPlot.AxisY[0].MaxValue * 1500;

        }
        
        private void getPlotDataFromDataTree(string sequence, string protein, string file)
        {
            while (elutionPlot.Series.Count > 0)
                elutionPlot.Series.RemoveAt(elutionPlot.Series.Count() - 1);

            while (spectrumPlot.Series.Count > 0)
                spectrumPlot.Series.RemoveAt(spectrumPlot.Series.Count() - 1);

            while (spectrumPlot.VisualElements.Count() > 0)
                spectrumPlot.VisualElements.RemoveAt(spectrumPlot.VisualElements.Count() - 1);


            var rawFilePath = "";
            foreach (DataGridViewRow row in dataUpload.Rows)
            {
                if (Path.GetFileNameWithoutExtension(row.Cells[0].Value.ToString()).Equals(file))
                {
                    rawFilePath = row.Cells[1].Value.ToString();
                }
            }

            var sqlPath = string.Format("Data Source={0}\\MyDatabase.sqlite; Version=3;", outputPath.Text);
            SQLiteConnection sqlReader = new SQLiteConnection(sqlPath);
            sqlReader.Open();

            var query = string.Format("SELECT * FROM AllGlycoPSMs WHERE AllGlycoPSMs.Sequence='{0}' AND AllGlycoPSMs.ProteinFasta='{1}' AND AllGlycoPSMs.File='{2}' ORDER BY AllGlycoPSMs.Score DESC", sequence, protein, file);
            var command = new SQLiteCommand(query, sqlReader);
            var reader = command.ExecuteReader();

            //For annotating MS Spectrum
            var pepSequence = "";
            var glycans = new List<string>();
            var mods = new List<string>();
            var charge = 0;
            var spectrumNumber = 0;
            var RTs = new List<string>();
            var elutionIntensities = new List<string>();
            var glycoPSM = new GlycoPSM();

            while (reader.Read())
            {
                glycoPSM = SQLiteStringToGlycoPSM(reader);

                pepSequence = reader["Peptide"].ToString();
                glycans = reader["Glycans"].ToString().Split(';').ToList();
                mods = reader["ModsVar"].ToString().Split(';').ToList();
                charge = Int32.Parse(reader["Charge"].ToString());
                spectrumNumber = Int32.Parse(reader["ScanNum"].ToString());
                file = reader["File"].ToString();
                RTs = reader["RetentionTimes"].ToString().Split(',').ToList();
                elutionIntensities = reader["Intensities"].ToString().Split(',').ToList();

                break;
            }

            // Moved to separate method to neaten up
            elutionPlot.Series = ChangeElutionPlot(RTs, elutionIntensities, file, spectrumNumber);

            var fragFinder = new glycoFragmentFinder(rawFilePath, glycoPSM);
            fragFinder.ReturnData += HandleFragmentMatchDataReturn;
            Task thread = new Task(fragFinder.crunch);
            thread.Start();


        }
        
        private void HandleFragmentMatchDataReturn(object sender, FragmentDataReturnArgs e)
        {
            switch (e.returnType)
            {
                case "UpdateSpectrumPlot":
                    UpdateSpectrumPlotWithFragDataReturn(e);
                    break;

                case "BatchAnnotationUpdate":
                    AddProgressText(string.Format("Spectra from {0} annotated in {1}", e.file, e.timeTaken));
                    break;
            }

        }

        private void UpdateSpectrumPlotWithFragDataReturn(FragmentDataReturnArgs data)
        {
            if (InvokeRequired)
            {
                spectrumPlot.Invoke(new Action<FragmentDataReturnArgs>(UpdateSpectrumPlotWithFragDataReturn), data);
                return;
            }

            var intensities = data.spectrum.GetIntensities();
            var mzs = data.spectrum.GetMasses();

            spectrumPlot.Series = new SeriesCollection();
            spectrumPlot.VisualElements.Clear();

            //Live Charts Implementation
            //var data = new ChartValues<ObservablePoint>();
            var plotData = new ChartValues<ObservablePoint>();


            for (int i = 0; i < mzs.Count(); i++)
            {
                plotData.Add(new ObservablePoint
                {
                    X = mzs[i],
                    Y = intensities[i],

                });
            }

            spectrumPlot.Series.Add(new ColumnSeries
            {
                Values = plotData,
                ColumnPadding = -0.2,
                Stroke = System.Windows.Media.Brushes.Black,
                StrokeThickness = 1000,
                Width = 100,

            });


            // Return is [0] List of intact peptide frags [1] Neutral loss peptide frags [2] PepHexMatches [3] Y ions [4] oxonium ions
            //var spectrumData = GetSpectrumAnnotations(pepSequence, spectrum, mods, glycans, charge);

            var spectrumData = new List<List<FragmentMatch>>() { data.peptideFragments, data.peptideNeutralLossFragments, data.peptideFragmentsMustIncludeGlycan, data.YIons, data.oxoniumIons };

            //Used for testing ScottPlot
            //for(int i = 0; i < spectrumData.Count(); i++)
            //{
            //    var writer = new StreamWriter(string.Format(@"{0}\SpectrumPlot_{1}.csv", outputPath.Text, i));

            //    writer.WriteLine("MZ,Intensity,Label");
            //    foreach(var frag in spectrumData[i])
            //    {
            //        writer.WriteLine("{0},{1},{2}", frag.fragmentMZ, frag.fragmentSignal, frag.label);
            //    }

            //    writer.Close();

            //}

            //var writer2 = new StreamWriter(string.Format(@"{0}\SpectrumPlot.csv", outputPath.Text));
            //writer2.WriteLine("MZ,Intensity,Label");
            //for (int i = 0; i < intensities.Count(); i++)
            //{
            //    writer2.WriteLine("{0},{1},", mzs[i], intensities[i]);
            //}
            //writer2.Close();

            var brushes = new List<System.Windows.Media.Brush>() { System.Windows.Media.Brushes.Red,
                                                                   System.Windows.Media.Brushes.Blue,
                                                                   System.Windows.Media.Brushes.Green,
                                                                   System.Windows.Media.Brushes.Brown,
                                                                   System.Windows.Media.Brushes.DarkGray};
            var j = 0;
            foreach (var spec in spectrumData)
            {

                var specIntensities = spec.Select(x => x.fragmentSignal).ToList();
                var specMZs = spec.Select(x => x.fragmentMZ).ToList();
                var labels = spec.Select(x => x.label).ToList();

                var specData = new ChartValues<ObservablePoint>();
                for (int i = 0; i < specMZs.Count(); i++)
                {
                    specData.Add(new ObservablePoint
                    {
                        X = specMZs[i],
                        Y = specIntensities[i]
                    });

                    var label = new VisualElement
                    {
                        X = specMZs[i] + 10,
                        Y = specIntensities[i],
                        VerticalAlignment = System.Windows.VerticalAlignment.Top,
                        HorizontalAlignment = System.Windows.HorizontalAlignment.Right,
                        UIElement = new TextBlock { Text = labels[i], RenderTransformOrigin = new System.Windows.Point(0.0, 0.9), RenderTransform = new RotateTransform(-90) }, // , LayoutTransform = new RotateTransform(-90)},                        
                        //UIElement.RenderTransform = new RotateTransform(-90),
                        //LayoutTransform = new RotateTransform(-90),
                    };
                    spectrumPlot.VisualElements.Add(label);
                }

                spectrumPlot.Series.Add(new ColumnSeries
                {
                    Values = specData,
                    ColumnPadding = -2,
                    Stroke = brushes[j],
                    StrokeThickness = 1000,
                    Width = 1000,

                });

                j++;
            }

            //Zed Graph Implementation
            /**
            GraphPane barPane = spectrumPlot.GraphPane;
            BarItem curve = barPane.AddBar("Spectrum", mzs, intensities, Color.Black);
            spectrumPlot.AxisChange();
            spectrumPlot.Refresh()
            **/

        }

        private void BatchAnnotateSpectra()
        {
            var sqlPath = string.Format("Data Source={0}\\MyDatabase.sqlite; Version=3;", outputPath.Text);
            SQLiteConnection sqlReader = new SQLiteConnection(sqlPath);
            sqlReader.Open();

            var uniqueFiles = new List<string>();

            var fileQuery = "SELECT DISTINCT File FROM AllGlycoPSMS";
            var fileCommand = new SQLiteCommand(fileQuery, sqlReader);
            var fileReader = fileCommand.ExecuteReader();

            while (fileReader.Read())
            {
                uniqueFiles.Add(fileReader["File"].ToString());
            }

            //var tasks = new List<Task>();            
            foreach (var file in uniqueFiles)
            {
                AddProgressText(string.Format("Annotating spectra from {0}", file));
                var rawFilePath = "";
                foreach (DataGridViewRow row in dataUpload.Rows)
                {
                    if (Path.GetFileNameWithoutExtension(row.Cells[0].Value.ToString()).Equals(file))
                    {
                        rawFilePath = row.Cells[1].Value.ToString();
                    }
                }

                var query = string.Format("SELECT * FROM AllGlycoPSMs WHERE AllGlycoPSMs.File='{0}'", file);
                var command = new SQLiteCommand(query, sqlReader);
                var reader = command.ExecuteReader();

                var psms = new List<GlycoPSM>();
                while (reader.Read())
                {
                    var glycoPSM = SQLiteStringToGlycoPSM(reader);
                    psms.Add(glycoPSM);
                }

                //var watch = new Stopwatch();
                //watch.Start();
                var annotator = new glycoFragmentFinder(rawFilePath, psms, outputPath.Text, PerformGlycanLocalization.Checked);
                annotator.ReturnData += HandleFragmentMatchDataReturn;
                annotator.useType = "BatchAnnotation";
                //var task = new Task(annotator.crunch);
                //Task.WaitAll(task);
                //var task = new Task(annotator.crunch);
                //tasks.Add(task);
                //task.Start();
                annotator.crunch();
                //watch.Stop();
                //AddProgressText(string.Format("Spectra annotated in {0}", watch.Elapsed));

                //task.Start();
            }
            //Task.WaitAll(tasks.ToArray());
        }
        
        private void dataTree_AfterSelect(object sender, TreeViewEventArgs e)
        {
            if (e.Node.Level == 2)
            {
                var sequence = e.Node.Text;
                var prot = e.Node.Parent.Text;
                var file = e.Node.Parent.Parent.Text;
                getPlotDataFromDataTree(sequence, prot, file);
            }
        }    

        // Return is [0] List of intact peptide frags [1] Neutral loss peptide frags [2] PepHexMatches [3] Y ions [4] oxonium ions
        private List<List<SpectrumFragment>> GetSpectrumAnnotations(string pep, ThermoSpectrum spectrum, List<string> mods, List<string> glycans, int charge)
        {
            var peptide = new Peptide(pep);
            var peptideNoGlycans = new Peptide(pep);
            var peptideHexNAc = new Peptide(pep);

            var hexnac = new ChemicalFormula("C8O5NH13");

            foreach (var mod in mods)
            {
                // e.g. M3(Oxidation / 15.9949) and N5(NGlycan / 1378.4756)
                var modMass = Double.Parse(mod.Split(' ')[2].Split(')')[0]);
                var modPosition = Int32.Parse(mod.Split('(')[0].Substring(1));
                var modName = mod.Split('(')[1].Split(' ')[0];

                var newMod = new Modification(modMass, modName);

                if (mod.Contains("Glycan"))
                {
                    peptide.AddModification(newMod, modPosition);
                    peptideHexNAc.AddModification(new Modification(hexnac.MonoisotopicMass, "HexNAc"), modPosition);
                }
                else
                {
                    peptide.AddModification(newMod, modPosition);
                    peptideHexNAc.AddModification(newMod, modPosition);
                    peptideNoGlycans.AddModification(newMod, modPosition);
                }

            }

            var peptideFragments = GeneratePeptideFragments(peptide);
            var peptideMatches = LookForPeptideFragments(peptideFragments, spectrum, charge);

            var peptideNoGlycanFragments = GeneratePeptideFragments(peptideNoGlycans);
            var peptideNoGlycanMatches = LookForPeptideFragments(peptideNoGlycanFragments, spectrum, charge);

            var peptideHexNAcFragments = GeneratePeptideFragments(peptideHexNAc);
            var peptideHexNAcMatches = LookForPeptideFragments(peptideHexNAcFragments, spectrum, charge);

            var YIonMatches = LookForYIonFragments(peptideNoGlycans, glycans, spectrum, charge);

            var oxoniumIonMathces = LookForOxoniumIons(glycans, spectrum);

            var returnList = new List<List<SpectrumFragment>>() { peptideMatches, peptideNoGlycanMatches, peptideHexNAcMatches, YIonMatches, oxoniumIonMathces };

            return returnList;

        }

        private List<SpectrumFragment> LookForOxoniumIons(List<string> glycans, ThermoSpectrum spectrum)
        {

            var returnList = new List<SpectrumFragment>();

            foreach (var gly in glycans)
            {
                var glycan = new Glycan(gly);
                glycan.GenerateCombinations();
                var glycanPieces = glycan.AllFragments;
                var filteredPieces = GlycanMethods.GetValidGlycanStructures(glycanPieces);

                foreach (var piece in filteredPieces)
                {
                    var matches = new List<ThermoMzPeak>();
                    for (int i = 0; i < 4; i++)
                    {
                        var mz = GlycanMethods.GetGlycanMass(piece.Value) + (i * Constants.Hydrogen);

                        var outpeaks = new List<ThermoMzPeak>();
                        if (spectrum.TryGetPeaks(DoubleRange.FromPPM(mz, 40), out outpeaks))
                        {
                            matches.AddRange(outpeaks);
                        }
                    }
                    if (matches.Count() > 0)
                    {
                        var mostIntense = matches.OrderByDescending(x => x.Intensity).ToList()[0];
                        returnList.Add(new SpectrumFragment(mostIntense.MZ, mostIntense.Intensity, piece.Key));
                    }
                }
            }

            return returnList;
        }

        private List<SpectrumFragment> LookForYIonFragments(Peptide pep, List<string> glycans, ThermoSpectrum spectrum, int charge)
        {
            var returnList = new List<SpectrumFragment>();

            foreach (var gly in glycans)
            {
                var glycan = new Glycan(gly);
                glycan.GenerateCombinations();
                var glycanPieces = glycan.AllFragments;
                var filteredPieces = GlycanMethods.GetValidGlycanStructures(glycanPieces);

                foreach (var piece in filteredPieces)
                {
                    for (int i = charge; i > 0; i--)
                    {
                        var matches = new List<ThermoMzPeak>();
                        double YIonMZ = ((pep.MonoisotopicMass) + (Constants.Hydrogen * ((double)i)) + GlycanMethods.GetGlycanMass(piece.Value)) / ((double)i);
                        if (YIonMZ < 3000)
                        {
                            for (int j = 0; j < 4; j++)
                            {
                                var isotopeMZ = (YIonMZ + (j * Constants.Hydrogen)) / ((double)i);
                                var outpeaks = new List<ThermoMzPeak>();
                                if (spectrum.TryGetPeaks(DoubleRange.FromPPM(isotopeMZ, 40), out outpeaks))
                                {
                                    matches.AddRange(outpeaks);
                                }
                            }
                        }
                        if (matches.Count() > 0)
                        {
                            var mostIntense = matches.OrderByDescending(x => x.Intensity).ToList()[0];
                            returnList.Add(new SpectrumFragment(mostIntense.MZ, mostIntense.Intensity, "Pep+" + piece.Key + "_z" + i.ToString()));
                        }
                    }
                }
            }

            return returnList;
        }

        private List<SpectrumFragment> LookForPeptideFragments(List<Fragment> frags, ThermoSpectrum spectrum, int charge)
        {
            var returnPeaks = new List<SpectrumFragment>();

            foreach (var frag in frags)
            {
                for (int i = 1; i < charge; i++)
                {
                    var fragMZ = (frag.MonoisotopicMass + (Constants.Hydrogen * (double)i)) / (double)i;

                    var outpeaks = new List<ThermoMzPeak>();
                    if (spectrum.TryGetPeaks(DoubleRange.FromPPM(fragMZ, 15), out outpeaks))
                    {
                        var mostIntense = outpeaks.OrderBy(x => x.Intensity).ToList()[0];

                        returnPeaks.Add(new SpectrumFragment(mostIntense.MZ, mostIntense.Intensity, frag.Type.ToString() + frag.Number.ToString() + "_z" + i.ToString()));
                    }
                }
            }

            return returnPeaks;
        }

        private List<Fragment> GeneratePeptideFragments(Peptide peptide)
        {
            var pepFrags = peptide.Fragment(FragmentTypes.b).ToList().ToList();
            pepFrags.AddRange(peptide.Fragment(FragmentTypes.c));
            pepFrags.AddRange(peptide.Fragment(FragmentTypes.y));
            pepFrags.AddRange(peptide.Fragment(FragmentTypes.zdot));

            return pepFrags;
        }
        
        private void dataTree_MouseUp(object sender, MouseEventArgs e)
        {

        }

        private void dataTree_MouseDown(object sender, MouseEventArgs e)
        {

        }

        private void dataTree_MouseMove(object sender, MouseEventArgs e)
        {

        }

        private void splitContainer5_Panel1_Paint(object sender, PaintEventArgs e)
        {

        }

        private void splitContainer4_Panel2_Paint(object sender, PaintEventArgs e)
        {

        }

        private void splitContainer3_Panel1_Resize(object sender, EventArgs e)
        {

        }

        private void elutionPlot_DataClick(object sender, ChartPoint chartPoint)
        {

        }

        private void processPSMs()
        {

        }    

        //Adds Row numbers to DataGridView
        private void dataGridView1_RowPostPaint(object sender, DataGridViewRowPostPaintEventArgs e)
        {
            var grid = sender as DataGridView;
            var rowIdx = (e.RowIndex + 1).ToString();

            var centerFormat = new StringFormat()
            {
                Alignment = StringAlignment.Center,
                LineAlignment = StringAlignment.Center
            };
            
            var headerBounds = new Rectangle(e.RowBounds.Left, e.RowBounds.Top, grid.RowHeadersWidth, e.RowBounds.Height);
            e.Graphics.DrawString(rowIdx, this.Font, SystemBrushes.ControlText, headerBounds, centerFormat);
        }

        private GlycoPSM SQLiteStringToGlycoPSM(SQLiteDataReader reader)
        {
            var sequence = reader["Sequence"].ToString();
            var peptidesToBeParsed = reader["PeptideParseFriendly"].ToString();
            var peptideStartPosition = Int32.Parse(reader["Position"].ToString());
            var modsToBeParsed = reader["ModsVar"].ToString();
            var glycansToBeParsed = reader["Glycans"].ToString();
            var PEP2D = double.Parse(reader["PEP2D"].ToString());
            var logProb = double.Parse(reader["logProb"].ToString());
            var score = double.Parse(reader["Score"].ToString());
            var deltaScore = double.Parse(reader["DeltaScore"].ToString());
            var deltaModScore = double.Parse(reader["DeltaModScore"].ToString());
            var charge = Int32.Parse(reader["Charge"].ToString());
            var mzObs = double.Parse(reader["ObsMZ"].ToString());
            var ppmError = double.Parse(reader["ppmerr"].ToString());
            var obsMH = double.Parse(reader["ObsMH"].ToString());
            var cleavage = reader["Cleavage"].ToString();
            var glycanPositions = reader["GlycansPos"].ToString();
            var proteinName = reader["ProteinFasta"].ToString();
            var scanTime = double.Parse(reader["ScanTime"].ToString());
            var scanNumber = Int32.Parse(reader["ScanNum"].ToString());
            var modsFixed = reader["ModsFixed"].ToString();
            var FDR2D = double.Parse(reader["FDR2D"].ToString());
            var FDR2Dunique = double.Parse(reader["FDRuniq2d"].ToString());
            var qvalue2D = double.Parse(reader["qValue2D"].ToString());
            var isGlycoPeptide = true;
            //var test1 = reader["isGlycoPeptide"].ToString(); // 'True' in SQLite Table reads in '0', should be '1'. Not a big issue because all psms are already filtered for glycopeptides at this point
            //var test = reader["isGlycoPeptide"]; //.Equals("True") ? true : false;
            var fragmentation = reader["DissociationType"].ToString();

            List<int> glycanPositionsList = new List<int>();
            List<string> mods = new List<string>();
            List<string> glycans = new List<string>();

            bool seenWithHCD = false;
            bool seenWithETD = false;
            bool NXSmotif = false;
            bool NXTmotif = false;
            bool isLocalized = false;
            bool Nlinked = false;
            bool Olinked = false;

            List<double> glycanMasses = new List<double>();

            bool matchedToUniprot = false;
            string uniprotEvidenceType = "None";
            int uniprotEvidenceNumber = 0;

            string[] parsedPeptide = peptidesToBeParsed.Split(',');
            Peptide peptide = new Peptide(parsedPeptide[0]);
            Peptide peptideNoGlycan = new Peptide(parsedPeptide[0]);                
            double peptideMonoMass = peptide.MonoisotopicMass;

            char[] peptideTermini = parsedPeptide[1].ToCharArray();
            char peptideCterminusNextResidue = peptideTermini[1];

            string[] parsedSequenceWithMods = sequence.Split('.');
            string sequenceWithMods = parsedSequenceWithMods[1];

            string[] parsedProteinName = proteinName.Split('|');
            string uniprotID = parsedProteinName[1];

            if (fragmentation.Equals("HCD"))
                seenWithHCD = true;

            if (fragmentation.Equals("ETD"))
                seenWithETD = true;

            if (!String.IsNullOrEmpty(glycanPositions))
            {
                string[] glycansPosParsedArray = glycanPositions.Split(';');
                for (int i = 0; i < glycansPosParsedArray.Length; i++)
                {
                    int glycanPos = Convert.ToInt32(glycansPosParsedArray[i]);
                    glycanPositionsList.Add(glycanPos);
                }
            }

            if (!String.IsNullOrEmpty(glycansToBeParsed))
            {
                string[] glycansParsedArrary = glycansToBeParsed.Split(';');

                for (int i = 0; i < glycansParsedArrary.Length; i++)
                {
                    string glycan = glycansParsedArrary[i];
                    if (glycan[0].Equals(' '))
                    {
                        glycan = glycan.Substring(1);
                    }
                    glycans.Add(glycan);
                }
            }

            int numberOfSites = glycans.Count;

            if (!String.IsNullOrEmpty(modsToBeParsed))
            {
                string[] modsParsedArrary = modsToBeParsed.Split(';');
                for (int i = 0; i < modsParsedArrary.Length; i++)
                {
                    string mod = modsParsedArrary[i];
                    if (mod[0].Equals(' '))
                    {
                        mod = mod.Substring(1);
                    }

                    mods.Add(mod);
                    string modName = GetModName(mod);
                    double modMass = GetModMass(mod);
                    int modPosition = GetModPosition(mod);
                    if (modName.Contains("Glycan"))
                    {
                        if (!peptide.Sequence[modPosition - 1].Equals('N'))
                        {
                            modName = "OGlycan";
                        }
                        else
                        {
                            modName = "NGlycan";
                        }
                        modName = modName + "_" + modMass;
                        glycanMasses.Add(modMass);
                    }

                    if (modName.Contains("NGlycan"))
                        Nlinked = true;

                    if (modName.Contains("OGlycan"))
                        Olinked = true;

                    Modification modToAdd = new Modification(modMass, modName);
                    peptide.AddModification(modToAdd, modPosition);

                    if (modName.Contains("NGlycan"))
                    {
                        if ((modPosition + 2) > peptide.Length)
                        {
                            if (!peptide.GetResidue(peptide.Length - 1).Equals('P'))
                            {
                                if (peptideCterminusNextResidue.Equals('S'))
                                {
                                    NXSmotif = true;
                                }
                                if (peptideCterminusNextResidue.Equals('T'))
                                {
                                    NXTmotif = true;
                                }
                            }
                        }
                        else
                        {
                            if (!peptide.GetResidue(modPosition).Letter.Equals('P'))
                            {
                                if (peptide.GetResidue(modPosition + 1).Letter.Equals('S'))
                                {
                                    NXSmotif = true;
                                }
                                if (peptide.GetResidue(modPosition + 1).Letter.Equals('T'))
                                {
                                    NXTmotif = true;
                                }
                            }
                        }
                    }
                }
            }

            if (!String.IsNullOrEmpty(modsFixed))
            {
                string[] modsParsedArrary = modsFixed.Split(';');
                for (int i = 0; i < modsParsedArrary.Length; i++)
                {
                    string mod = modsParsedArrary[i];
                    if (mod[0].Equals(' '))
                    {
                        mod = mod.Substring(1);
                    }
                    mods.Add(mod);
                    string modName = GetModName(mod);
                    double modMass = GetModMass(mod);
                    int modPosition = GetModPosition(mod);

                    Modification modToAdd = new Modification(modMass, modName);
                    peptide.AddModification(modToAdd, modPosition);
                }
            }

            if (deltaModScore >= 10)
                isLocalized = true;

            GlycoPSM psm = new GlycoPSM(peptide, peptideMonoMass, peptide.SequenceWithModifications, mods, glycans, glycanMasses, glycanPositionsList,
                uniprotID, PEP2D, logProb, score, deltaScore, deltaModScore, mzObs, charge, numberOfSites, ppmError, obsMH, cleavage, proteinName,
                peptideStartPosition, scanTime, scanNumber, FDR2D, FDR2Dunique, qvalue2D, fragmentation, isGlycoPeptide, seenWithHCD, seenWithETD,
                NXSmotif, NXTmotif, isLocalized, Nlinked, Olinked, matchedToUniprot, uniprotEvidenceType, uniprotEvidenceNumber);

            psm.file = reader["File"].ToString();
            psm.Condition = reader["Condition"].ToString();
            psm.Replicate = reader["Replicate"].ToString();

            return psm;
            
        }

        public static double GetModMass(string mod)
        {
            string[] modParse1 = mod.Split('(');
            string[] modParse2 = modParse1[1].Split(')');
            string[] modParse3 = modParse2[0].Split('/');

            string massSubstring = modParse3[1].Substring(1);

            double modMass = Convert.ToDouble(massSubstring);

            return modMass;
        }

        public static int GetModPosition(string mod)
        {
            string[] modParse1 = mod.Split('(');

            string modPositionString = modParse1[0].Substring(1);

            int modPosition = Convert.ToInt32(modPositionString);

            return modPosition;
        }

        public static string GetModName(string mod)
        {
            string[] modParse1 = mod.Split('(');
            string[] modParse2 = modParse1[1].Split(')');
            string[] modParse3 = modParse2[0].Split('/');
            string[] modParse4 = modParse3[0].Split(' ');

            string modName = modParse4[0];

            return modName;
        }

        private static SeriesCollection ChangeElutionPlot(List<string> RTs, List<string> elutionIntensities, string file, int scan)
        {
            var lineSeries = new SeriesCollection();
            var elutionData = new ChartValues<ObservablePoint>();
            for (int i = 0; i < RTs.Count() - 1; i++)
            {
                elutionData.Add(new ObservablePoint
                {
                    X = Double.Parse(RTs[i]),
                    Y = Double.Parse(elutionIntensities[i])

                });
            }

            lineSeries.Add(new LineSeries
            {
                Title = Path.GetFileNameWithoutExtension(file) + "_" + scan,
                Values = elutionData,
                PointGeometry = DefaultGeometries.Circle,
                PointGeometrySize = 5
            });
            return lineSeries;
        }

        private void PopulateGlycanPosComboBox(string peptideSequence, string originalPositionString)
        {
            glycanPosComboBox.Items.Clear();
            var position = Int32.Parse(originalPositionString.Split(';')[0]);
            glycanPosComboBox.Items.Add(peptideSequence[position - 1] + position.ToString() + "*");
            glycanPosComboBox.SelectedIndex = 0;

            for(int i = 0; i < peptideSequence.Length; i++)
            {
                try
                {
                    if ("ST".Contains(peptideSequence[i]) && (i + 1) != position)
                    {
                        glycanPosComboBox.Items.Add(peptideSequence[i] + (i + 1).ToString());
                    }

                    if ("N".Contains(peptideSequence[i]) && (i + 1) != position && "ST".Contains(peptideSequence[i + 2]) && !peptideSequence[i + 1].Equals("P"))
                    {
                        glycanPosComboBox.Items.Add(peptideSequence[i] + (i + 1).ToString());
                    }
                }catch(Exception e)
                {

                }
                
            }
        }

        private void viewPeptidesDataGrid_SelectionChanged(object sender, EventArgs e)
        {
            if(viewPeptidesDataGrid.SelectedRows.Count == 1 && dataUploaded)
            {
                var scanNumIndex = -1;
                var fileIndex = -1;
                var sequenceIndex = -1;
                var modsIndex = -1;
                
                foreach(DataGridViewColumn column in viewPeptidesDataGrid.Columns)
                {
                    if (column.HeaderText.Equals("ScanNum"))
                    {
                        scanNumIndex = column.Index;
                    }

                    if (column.HeaderText.Equals("File"))
                    {
                        fileIndex = column.Index;
                    }

                    if (column.HeaderText.Equals("Peptide"))
                    {
                        sequenceIndex = column.Index;
                    }

                    if (column.HeaderText.Equals("ModsVar"))
                    {
                        modsIndex = column.Index;
                    }
                    
                }                

                if(scanNumIndex != -1 && fileIndex != -1)
                {
                    var row = viewPeptidesDataGrid.SelectedRows[0];
                    var scanNum = Int32.Parse(row.Cells[scanNumIndex].Value.ToString());
                    var file = row.Cells[fileIndex].Value.ToString();
                    getPlotDataFromSQL(scanNum, file);
                }

                if(sequenceIndex != -1 && modsIndex != -1)
                {
                    var row = viewPeptidesDataGrid.SelectedRows[0];
                    var mods = row.Cells[modsIndex].Value.ToString().Split(';');
                    var position = 0;
                    foreach(var mod in mods)
                    {
                        if (mod.Contains("Glycan"))
                        {
                            position = Int32.Parse(mod.Split('(')[0].Substring(1));
                        }
                    }
                    PopulateGlycanPosComboBox(row.Cells[sequenceIndex].Value.ToString(), position.ToString());
                }                    
            }
        }

        //First list binMid, Second list binCount
        private List<double[]> GetBinCountsFromSqlReader(int numBins, SQLiteDataReader reader, double min, double max, string key)
        {            
            double binWidth = (max - min) / Convert.ToDouble(numBins);
                
            var binMids = new double[numBins];
            var binCounts = new double[numBins];            

            for(int i = 0; i < numBins; i++)
            {
                var binMin = (binWidth * Convert.ToDouble(i)) + min;
                var binMax = (binWidth * Convert.ToDouble(i + 1)) + min;
                binMids[i] = binMin + (binWidth / 2.0);
            }

            while (reader.Read())
            {
                var value = double.Parse(reader[key].ToString());

                for (int i = 0; i < numBins; i++)
                {
                    var binMin = (binWidth * Convert.ToDouble(i)) + min;
                    var binMax = (binWidth * Convert.ToDouble(i + 1)) + min;

                    if (value > binMin && value <= binMax)
                        binCounts[i]++;
                    
                }
            }           

            return new List<double[]>() { binMids, binCounts };
        }

        private void statsComboBox_SelectedValueChanged(object sender, EventArgs e)
        {
            if (dataUploaded)
            {
                var selection = statsComboBox.SelectedItem.ToString();
                if(selection.Equals("All Files"))
                {
                    fillFragHistPlot("", true);
                    fillSeqCoveragePlot("", true);
                    fillPrecursorMassErrorPlot("", true);
                    fillGlycanTypePieChart("", true);
                }
                else
                {
                    fillFragHistPlot(statsComboBox.SelectedItem.ToString(), false);
                    fillSeqCoveragePlot(statsComboBox.SelectedItem.ToString(), false);
                    fillPrecursorMassErrorPlot(statsComboBox.SelectedItem.ToString(), false);
                    fillGlycanTypePieChart(statsComboBox.SelectedItem.ToString(), false);
                }                
            }
        }

        private void viewPeptidesDataGrid_RowPostPaint(object sender, DataGridViewRowPostPaintEventArgs e)
        {
            var grid = sender as DataGridView;
            var rowIdx = (e.RowIndex + 1).ToString();

            var centerFormat = new StringFormat()
            {
                Alignment = StringAlignment.Center,
                LineAlignment = StringAlignment.Center
            };

            var headerBounds = new Rectangle(e.RowBounds.Left, e.RowBounds.Top, grid.RowHeadersWidth, e.RowBounds.Height);
            e.Graphics.DrawString(rowIdx, this.Font, SystemBrushes.ControlText, headerBounds, centerFormat);
        }

        private void resultViewBox_SelectedValueChanged(object sender, EventArgs e)
        {
            if (dataUploaded)
            {
                var condition = conditionComboBox.SelectedItem.ToString();
                var replicate = replicateComboBox.SelectedItem.ToString();
                var table = resultViewBox.SelectedItem.ToString();

                if(condition.Equals("All") && replicate.Equals("All"))
                {
                    fillResultTable(string.Format("SELECT * FROM All{0}", resultViewBox.SelectedItem.ToString()));
                    return;
                }                    

                if(!condition.Equals("All") && !replicate.Equals("All"))
                {
                    fillResultTable(string.Format("SELECT * FROM {0}_Condition_{1}_Replicate_{2}", table, condition, replicate));
                    return;
                }

                if(!condition.Equals("All") && replicate.Equals("All"))
                {
                    fillResultTable(string.Format("SELECT * FROM {0}_Condition_{1}", table, condition));
                }
            }
        }

        private void conditionComboBox_SelectedIndexChanged(object sender, EventArgs e)
        {
            var conditionText = conditionComboBox.SelectedItem.ToString();
            var tableText = resultViewBox.SelectedItem.ToString();

            replicateComboBox.Items.Clear();
            replicateComboBox.ResetText();           
            replicateComboBox.Items.Add("All");
            replicateComboBox.SelectedIndex = 0;

            if (conditionText.Equals("All"))
            {
                fillResultTable(string.Format("SELECT * FROM All{0}", tableText));
            }
            else
            {
                var connection = new SQLiteConnection("Data Source=" + outputPath.Text + "\\MyDatabase.sqlite; Version=3;");
                connection.Open();

                var command = new SQLiteCommand(string.Format("SELECT DISTINCT Replicate FROM AllGlycoPSMs WHERE AllGlycoPSMs.Condition='{0}'", conditionText), connection);
                var reader = command.ExecuteReader();

                while (reader.Read())
                {
                    var table = reader["Replicate"].ToString();
                    replicateComboBox.Items.Add(table);
                }

                connection.Close();
                GC.Collect();

                fillResultTable(string.Format("SELECT * FROM {0}_Condition_{1}", tableText, conditionText));
            }
        }

        private void replicateComboBox_SelectedIndexChanged(object sender, EventArgs e)
        {
            var replicateText = replicateComboBox.SelectedItem.ToString();
            var conditionText = conditionComboBox.SelectedItem.ToString();
            var tableText = resultViewBox.SelectedItem.ToString();
            
            if(replicateText.Equals("All") && conditionText.Equals("All"))
            {
                fillResultTable(string.Format("SELECT * FROM All{0}", tableText));
                return;
            }

            if(replicateText.Equals("All") && !conditionText.Equals("All"))
            {
                fillResultTable(string.Format("SELECT * FROM {0}_Condition_{1}", tableText, conditionText));
                return;
            }

            if(!replicateText.Equals("All") && !conditionText.Equals("All"))
            {
                fillResultTable(string.Format("SELECT * FROM {0}_Condition_{1}_Replicate_{2}", tableText, conditionText, replicateText));
            }
        }

        private void PerformProteinInference_CheckedChanged(object sender, EventArgs e)
        {
            if (PerformProteinInference.Checked)
            {
                PeptidesLabel.Visible = true;
                AdditionalPeptidesPath.Visible = true;
                ProteaseLabel.Visible = true;
                fastaLabel.Visible = true;
                fastaPath.Visible = true;
                ProteaseCB.Visible = true;
                miscleaveLabel.Visible = true;
                miscleaveCB.Visible = true;
                List<Protease> proteases = CSMSL.Proteomics.Protease.GetAllProteases().ToList();
                ProteaseCB.DataSource = proteases.Select(x => x.Name).ToList();
                var cleavages = new List<int> { 1, 2, 3, 4, 5 };
                miscleaveCB.DataSource = cleavages;
            }else
            {
                PeptidesLabel.Visible = false;
                AdditionalPeptidesPath.Visible = false;
                ProteaseLabel.Visible = false;
                fastaLabel.Visible = false;
                fastaPath.Visible = false;
                ProteaseCB.Visible = false;
                miscleaveCB.Visible = false;
                miscleaveLabel.Visible = false;
                ProteaseCB.DataSource = null;
                miscleaveCB.DataSource = null;
            }
        }

        private void AdditionalPeptidesPath_DragDrop(object sender, DragEventArgs e)
        {
            if (!e.Data.GetDataPresent(DataFormats.FileDrop))
                return;

            string[] files = e.Data.GetData(DataFormats.FileDrop) as string[];
            if (files == null)
                return;

            foreach (string file in files.Where(f => Path.GetExtension(f).Equals(".csv")))
            {
                AdditionalPeptidesPath.Text = file;
            }
        }

        private void AdditionalPeptidesPath_DragEnter(object sender, DragEventArgs e)
        {
            e.Effect = DragDropEffects.None;

            if (!e.Data.GetDataPresent(DataFormats.FileDrop))
                return;

            string[] files = e.Data.GetData(DataFormats.FileDrop) as string[];
            if (files == null)
                return;

            if (files.Any(f => Path.GetExtension(f).Equals(".csv")))
                e.Effect = DragDropEffects.Link;
        }

        private void fastaPath_DragDrop(object sender, DragEventArgs e)
        {
            if (!e.Data.GetDataPresent(DataFormats.FileDrop))
                return;

            string[] files = e.Data.GetData(DataFormats.FileDrop) as string[];
            if (files == null)
                return;

            foreach (string file in files.Where(f => Path.GetExtension(f).ToLower().Equals(".fasta")))
            {
                fastaPath.Text = file;
            }

            foreach (string file in files.Where(f => Path.GetExtension(f).ToLower().Equals(".fa")))
            {
                fastaPath.Text = file;
            }
        }

        private void fastaPath_DragEnter(object sender, DragEventArgs e)
        {
            e.Effect = DragDropEffects.None;

            if (!e.Data.GetDataPresent(DataFormats.FileDrop))
                return;

            string[] files = e.Data.GetData(DataFormats.FileDrop) as string[];
            if (files == null)
                return;

            if (files.Any(f => Path.GetExtension(f).ToLower().Equals(".fasta")))
                e.Effect = DragDropEffects.Link;

            if (files.Any(f => Path.GetExtension(f).ToLower().Equals(".fa")))
                e.Effect = DragDropEffects.Link;
        }

        private void label12_Click(object sender, EventArgs e)
        {

        }

        private void PerformGlycanLocalization_CheckedChanged(object sender, EventArgs e)
        {
            if (PerformGlycanLocalization.Checked)
            {
                var comboboxOptions = new int[] { 1 };
                this.modCountFilter.DataSource = comboboxOptions;
            }
            else
            {
                var comboboxOptions = new int[] { 1, 2, 3 };
                this.modCountFilter.DataSource = comboboxOptions;
            }
        }

        private void button2_Click(object sender, EventArgs e)
        {
            MessageBox.Show("This feature will analyze each PSM to test whether the MS2 spectra is better explained by changing the placement of the glycan modification. The feature is only available for singly-glycosylated peptides. Enabling glycan localization will increase the processing time by ~4-fold.");
        }

        private void button1_Click(object sender, EventArgs e)
        {
            MessageBox.Show("Byonic processes individual .raw files separately. Thus protein inference is not extended across an entire data set, resulting in cases of the same peptide sequence between result files having different protein assignments. This feature will perform protein inference across all Byonic result files, as wells as a list of additional peptide identifications provided as a .csv file. The additional peptides file must contain columns with 'Peptide' and 'PEP' headers. For example, one can provide peptide results from deamidated or whole proteome analyses. The additional peptides file is not required, however the user must provide a protein database in FASTA format and specify the correct protease and allowed number of miscleavages.");
        }

        private void button3_Click(object sender, EventArgs e)
        {
            MessageBox.Show("Loss of sugar monomers from a glycopeptide during ionization and ion transmission can result in artifacts of in-source fragmentation being identified as bonafine glycopeptides. This feature will determine if an identification is the result of an in-source fragment of a larger glycopeptide by searching raw files for co-elution of a larger species.");
        }

        private void inSourceFragData_SelectionChanged(object sender, EventArgs e)
        {
            if (inSourceFragData.SelectedRows.Count != 1 || !dataUploaded)
                return;

            var row = inSourceFragData.SelectedRows[0];

            var idRTIndex = -1;
            var idIntIndex = -1;
            var parentRTIndex = -1;
            var parentIntIndex = -1;

            foreach (DataGridViewColumn column in inSourceFragData.Columns)
            {

                if (column.HeaderText.Equals("ID_RTS"))
                    idRTIndex = column.Index;

                if (column.HeaderText.Equals("ID_INTENSITIES"))
                    idIntIndex = column.Index;

                if (column.HeaderText.Equals("PARENT_RTS"))
                    parentRTIndex = column.Index;

                if (column.HeaderText.Equals("PARENT_INTENSITIES"))
                    parentIntIndex = column.Index;

            }

            if(idRTIndex != -1 && idIntIndex != -1 && parentRTIndex != -1 && parentIntIndex != -1)
            {
                var idRTs = StringListToDoubleList(row.Cells[idRTIndex].Value.ToString().Split(';').ToList());
                var idInts = StringListToDoubleList(row.Cells[idIntIndex].Value.ToString().Split(';').ToList());
                var parentRTs = StringListToDoubleList(row.Cells[parentRTIndex].Value.ToString().Split(';').ToList());
                var parentInts = StringListToDoubleList(row.Cells[parentIntIndex].Value.ToString().Split(';').ToList());


                PlotInSourceData(idRTs, idInts, parentRTs, parentInts);

            }

        }

        private void PlotInSourceData(List<double> idRTs, List<double> idInts, List<double> parentRTs, List<double> parentInts)
        {
            if (InvokeRequired)
            {
                inSourceFragChart.Invoke(new Action<List<double>, List<double>, List<double>, List<double>>(PlotInSourceData), idRTs, idInts, parentRTs, parentInts);
                return;
            }

            while (inSourceFragChart.Series.Count > 0)
                inSourceFragChart.Series.RemoveAt(inSourceFragChart.Series.Count() - 1);


            var lineSeriesID = new SeriesCollection();
            var elutionDataID = new ChartValues<ObservablePoint>();

            var lineSeriesParent = new SeriesCollection();
            var elutionDataParent = new ChartValues<ObservablePoint>();

            for (int i = 0; i < idRTs.Count(); i++)
            {
                elutionDataID.Add(new ObservablePoint
                {
                    X = idRTs[i],
                    Y = idInts[i]                    
                });
            }

            for (int i = 0; i < parentRTs.Count(); i++)
            {
                elutionDataParent.Add(new ObservablePoint
                {
                    X = parentRTs[i],
                    Y = parentInts[i]
                });
            }

            lineSeriesID.Add(new LineSeries
            {
                Values = elutionDataID,
                PointGeometry = DefaultGeometries.Circle,
                PointGeometrySize = 5,
                LineSmoothness = 0.0               
            });

            lineSeriesID.Add(new LineSeries
            {
                Values = elutionDataParent,
                PointGeometry = DefaultGeometries.Circle,
                PointGeometrySize = 5,
                LineSmoothness = 0.0
            });

            inSourceFragChart.Series = lineSeriesID;
        }

        private List<double> StringListToDoubleList(List<string> list)
        {
            var returnList = new List<double>();

            foreach (var entry in list)
            {
                if (!string.IsNullOrEmpty(entry))
                    returnList.Add(double.Parse(entry));
            }

            return returnList;
        }

        private void inSourceFragData_RowPostPaint_1(object sender, DataGridViewRowPostPaintEventArgs e)
        {
            var grid = sender as DataGridView;
            var rowIdx = (e.RowIndex + 1).ToString();

            var centerFormat = new StringFormat()
            {
                Alignment = StringAlignment.Center,
                LineAlignment = StringAlignment.Center
            };

            var headerBounds = new Rectangle(e.RowBounds.Left, e.RowBounds.Top, grid.RowHeadersWidth, e.RowBounds.Height);
            e.Graphics.DrawString(rowIdx, this.Font, SystemBrushes.ControlText, headerBounds, centerFormat);
        }

        private void splitContainer7_Panel1_Paint(object sender, PaintEventArgs e)
        {

        }

        private GlycoPSM SQLiteStringToGlycoPSM(SQLiteDataReader reader, int glycanPosition)
        {
            var sequence = reader["Sequence"].ToString();
            var peptidesToBeParsed = reader["PeptideParseFriendly"].ToString();
            var peptideStartPosition = Int32.Parse(reader["Position"].ToString());
            var modsToBeParsed = reader["ModsVar"].ToString();
            var glycansToBeParsed = reader["Glycans"].ToString();
            var PEP2D = double.Parse(reader["PEP2D"].ToString());
            var logProb = double.Parse(reader["logProb"].ToString());
            var score = double.Parse(reader["Score"].ToString());
            var deltaScore = double.Parse(reader["DeltaScore"].ToString());
            var deltaModScore = double.Parse(reader["DeltaModScore"].ToString());
            var charge = Int32.Parse(reader["Charge"].ToString());
            var mzObs = double.Parse(reader["ObsMZ"].ToString());
            var ppmError = double.Parse(reader["ppmerr"].ToString());
            var obsMH = double.Parse(reader["ObsMH"].ToString());
            var cleavage = reader["Cleavage"].ToString();
            var glycanPositions = reader["GlycansPos"].ToString();
            var proteinName = reader["ProteinFasta"].ToString();
            var scanTime = double.Parse(reader["ScanTime"].ToString());
            var scanNumber = Int32.Parse(reader["ScanNum"].ToString());
            var modsFixed = reader["ModsFixed"].ToString();
            var FDR2D = double.Parse(reader["FDR2D"].ToString());
            var FDR2Dunique = double.Parse(reader["FDRuniq2d"].ToString());
            var qvalue2D = double.Parse(reader["qValue2D"].ToString());
            var isGlycoPeptide = true;
            //var test1 = reader["isGlycoPeptide"].ToString(); // 'True' in SQLite Table reads in '0', should be '1'. Not a big issue because all psms are already filtered for glycopeptides at this point
            //var test = reader["isGlycoPeptide"]; //.Equals("True") ? true : false;
            var fragmentation = reader["DissociationType"].ToString();

            List<int> glycanPositionsList = new List<int>();
            List<string> mods = new List<string>();
            List<string> glycans = new List<string>();

            bool seenWithHCD = false;
            bool seenWithETD = false;
            bool NXSmotif = false;
            bool NXTmotif = false;
            bool isLocalized = false;
            bool Nlinked = false;
            bool Olinked = false;

            List<double> glycanMasses = new List<double>();

            bool matchedToUniprot = false;
            string uniprotEvidenceType = "None";
            int uniprotEvidenceNumber = 0;

            string[] parsedPeptide = peptidesToBeParsed.Split(',');
            Peptide peptide = new Peptide(parsedPeptide[0]);
            Peptide peptideNoGlycan = new Peptide(parsedPeptide[0]);
            double peptideMonoMass = peptide.MonoisotopicMass;

            char[] peptideTermini = parsedPeptide[1].ToCharArray();
            char peptideCterminusNextResidue = peptideTermini[1];

            string[] parsedSequenceWithMods = sequence.Split('.');
            string sequenceWithMods = parsedSequenceWithMods[1];

            string[] parsedProteinName = proteinName.Split('|');
            string uniprotID = parsedProteinName[1];

            if (fragmentation.Equals("HCD"))
                seenWithHCD = true;

            if (fragmentation.Equals("ETD"))
                seenWithETD = true;

            if (!String.IsNullOrEmpty(glycanPositions))
            {
                string[] glycansPosParsedArray = glycanPositions.Split(';');
                for (int i = 0; i < glycansPosParsedArray.Length; i++)
                {
                    int glycanPos = Convert.ToInt32(glycansPosParsedArray[i]);
                    glycanPositionsList.Add(glycanPos);
                }
            }

            if (!String.IsNullOrEmpty(glycansToBeParsed))
            {
                string[] glycansParsedArrary = glycansToBeParsed.Split(';');

                for (int i = 0; i < glycansParsedArrary.Length; i++)
                {
                    string glycan = glycansParsedArrary[i];
                    if (glycan[0].Equals(' '))
                    {
                        glycan = glycan.Substring(1);
                    }
                    glycans.Add(glycan);
                }
            }

            int numberOfSites = glycans.Count;

            if (!String.IsNullOrEmpty(modsToBeParsed))
            {
                string[] modsParsedArrary = modsToBeParsed.Split(';');
                for (int i = 0; i < modsParsedArrary.Length; i++)
                {
                    string mod = modsParsedArrary[i];
                    if (mod[0].Equals(' '))
                    {
                        mod = mod.Substring(1);
                    }

                    mods.Add(mod);
                    string modName = GetModName(mod);
                    double modMass = GetModMass(mod);
                    
                    int modPosition = GetModPosition(mod);
                    if (modName.Contains("Glycan"))
                    {
                        // Hard code glycan mod position here
                        modPosition = glycanPosition;
                        if (!peptide.Sequence[modPosition - 1].Equals('N'))
                        {
                            modName = "OGlycan";
                        }
                        else
                        {
                            modName = "NGlycan";
                        }
                        modName = modName + "_" + modMass;
                        glycanMasses.Add(modMass);
                    }

                    if (modName.Contains("NGlycan"))
                        Nlinked = true;

                    if (modName.Contains("OGlycan"))
                        Olinked = true;

                    Modification modToAdd = new Modification(modMass, modName);
                    peptide.AddModification(modToAdd, modPosition);

                    if (modName.Contains("NGlycan"))
                    {
                        if ((modPosition + 2) > peptide.Length)
                        {
                            if (!peptide.GetResidue(peptide.Length - 1).Equals('P'))
                            {
                                if (peptideCterminusNextResidue.Equals('S'))
                                {
                                    NXSmotif = true;
                                }
                                if (peptideCterminusNextResidue.Equals('T'))
                                {
                                    NXTmotif = true;
                                }
                            }
                        }
                        else
                        {
                            if (!peptide.GetResidue(modPosition).Letter.Equals('P'))
                            {
                                if (peptide.GetResidue(modPosition + 1).Letter.Equals('S'))
                                {
                                    NXSmotif = true;
                                }
                                if (peptide.GetResidue(modPosition + 1).Letter.Equals('T'))
                                {
                                    NXTmotif = true;
                                }
                            }
                        }
                    }
                }
            }

            if (!String.IsNullOrEmpty(modsFixed))
            {
                string[] modsParsedArrary = modsFixed.Split(';');
                for (int i = 0; i < modsParsedArrary.Length; i++)
                {
                    string mod = modsParsedArrary[i];
                    if (mod[0].Equals(' '))
                    {
                        mod = mod.Substring(1);
                    }
                    mods.Add(mod);
                    string modName = GetModName(mod);
                    double modMass = GetModMass(mod);
                    int modPosition = GetModPosition(mod);

                    Modification modToAdd = new Modification(modMass, modName);
                    peptide.AddModification(modToAdd, modPosition);
                }
            }

            if (deltaModScore >= 10)
                isLocalized = true;

            GlycoPSM psm = new GlycoPSM(peptide, peptideMonoMass, peptide.SequenceWithModifications, mods, glycans, glycanMasses, glycanPositionsList,
                uniprotID, PEP2D, logProb, score, deltaScore, deltaModScore, mzObs, charge, numberOfSites, ppmError, obsMH, cleavage, proteinName,
                peptideStartPosition, scanTime, scanNumber, FDR2D, FDR2Dunique, qvalue2D, fragmentation, isGlycoPeptide, seenWithHCD, seenWithETD,
                NXSmotif, NXTmotif, isLocalized, Nlinked, Olinked, matchedToUniprot, uniprotEvidenceType, uniprotEvidenceNumber);

            psm.file = reader["File"].ToString();
            psm.Condition = reader["Condition"].ToString();
            psm.Replicate = reader["Replicate"].ToString();

            return psm;

        }

        private void glycanPosComboBox_SelectedIndexChanged(object sender, EventArgs e)
        {
            try
            {
                var text = glycanPosComboBox.Text;
                var position = Int32.Parse(text.Substring(1).Trim('*'));

                var row = viewPeptidesDataGrid.SelectedRows[0];

                var fileIndex = -1;
                var scanIndex = -1;

                foreach (DataGridViewColumn column in viewPeptidesDataGrid.Columns)
                {
                    if (column.HeaderText.Equals("ScanNum"))
                    {
                        scanIndex = column.Index;
                    }

                    if (column.HeaderText.Equals("File"))
                    {
                        fileIndex = column.Index;
                    }
                }

                if(fileIndex != -1 && scanIndex != -1)
                {
                    var connection = new SQLiteConnection(string.Format("Data Source={0}; Version=3;", outputPath.Text + "\\MyDatabase.sqlite"));
                    connection.Open();

                    var file = row.Cells[fileIndex].Value.ToString();
                    var scanNum = row.Cells[scanIndex].Value.ToString();

                    var rawFilePath = "";
                    foreach (DataGridViewRow row2 in dataUpload.Rows)
                    {
                        if (Path.GetFileNameWithoutExtension(row2.Cells[0].Value.ToString()).Equals(file))
                        {
                            rawFilePath = row2.Cells[1].Value.ToString();
                        }
                    }

                    var commandString = string.Format("SELECT * FROM AllGlycoPSMs WHERE AllGlycoPSMs.File='{0}' AND AllGlycoPSMs.ScanNum='{1}'", file, scanNum);
                    var command = new SQLiteCommand(commandString, connection);
                    var reader = command.ExecuteReader();

                    while (reader.Read())
                    {
                        var glyPsm = SQLiteStringToGlycoPSM(reader, position);
                        var fragFinder = new glycoFragmentFinder(rawFilePath, glyPsm);
                        fragFinder.useType = "SingleSpectrum";
                        fragFinder.ReturnData += HandleFragmentMatchDataReturn;
                        fragFinder.crunch();
                        peptideLabel.Text = glyPsm.peptide.GetSequenceWithModifications();
                        peptideLabel.Visible = true;
                    }

                    connection.Close();

                }
                

            }
            catch(Exception exc)
            {

            }            

        }
    }
}
