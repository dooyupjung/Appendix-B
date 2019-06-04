#define _CRT_SECURE_NO_WARNINGS

#define NAME 200
#define MAX 15000
#define LENGTH 3
#define AATYPE 26

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <direct.h>

int main(void) {
	FILE *fmotif_list, *fspace_list, *fmotif, *fspace, *fint, *fsum, *fname;
	char filename[NAME], seqname[NAME], newdir[NAME], whole[MAX], **motif, interaction[LENGTH * 2 + 1];
	int i, j, k, l, m, temp, len, motif_type, cnt, flag, total_int, total_type, num_proteins;

	fmotif_list = fopen("Raw_motifs/filename.txt", "r");
	fspace_list = fopen("Without_space/filename.txt", "r");

	system("mkdir Interaction_DB\\");
	fsum = fopen("Interaction_DB/summary.txt", "w");
	fname = fopen("Interaction_DB/filename.txt", "w");
	
	fprintf(fsum, "Summary of Interaction DB Construction:\n");

	while (1) {
		// Motif ���� ����Ʈ ��������
		for (i = 0; i < NAME; i++) {
			filename[i] = '\0';
			newdir[i] = '\0';
		}

		fgets(filename, NAME, fmotif_list);
		fscanf(fmotif_list, "%d\n", &num_proteins);

		for (i = 0; i < NAME; i++) {
			if (filename[i] == 10) {
				filename[i] = '\0';
				break;
			}
		}

		if (i == 200)	break;		// ���� ������ ����� �������� ����!

		strcat_s(newdir, NAME, "Raw_motifs/");
		strcat_s(newdir, NAME, filename);
		fmotif = fopen(newdir, "r");

		// ���� ���� Seq ���� ����Ʈ ��������
		for (i = 0; i < NAME; i++) {
			filename[i] = '\0';
			newdir[i] = '\0';
		}
		
		fgets(filename, NAME, fspace_list);
		fscanf(fspace_list, "%d\n", &num_proteins);
		filename[strlen(filename) - 1] = '\0';
		strcat_s(newdir, NAME, "Without_space/");
		strcat_s(newdir, NAME, filename);
		fspace = fopen(newdir, "r");
		
		// Interaction DB ���� ���� ����
		for (i = 0; i < NAME; i++)
			newdir[i] = '\0';
		strcat_s(newdir, NAME, "Interaction_DB/");
		filename[strlen(filename) - 5] = 'i';
		filename[strlen(filename) - 4] = 'n';
		filename[strlen(filename) - 3] = 't';
		filename[strlen(filename) - 2] = 'D';
		filename[strlen(filename) - 1] = 'B';
		fprintf(fname, "%s\n", filename);	// Interaction DB ���� ��� ������Ʈ!
		strcat_s(newdir, NAME, filename);
		fint = fopen(newdir, "w");
		
		fprintf(fsum, "\n");
		fprintf(fname, "%d\n", num_proteins);

		while (1) {
			total_int = 0;
			total_type = 0;

			// Seqname �Է� �ޱ�
			for (j = 0; j < NAME; j++)
				seqname[j] = '\0';
			fgets(seqname, NAME, fmotif);
			if (seqname[0] == '\0')	break;	// ���̻� seqname ������ ����!

			for (j = 0; j < NAME; j++)
				seqname[j] = '\0';
			fgets(seqname, NAME, fspace);
			fprintf(fsum, "%s", &seqname[1]);
			printf("%s", &seqname);

			// Seqname�� Interation_DB ���� ���� ���
			fprintf(fint, "%s", seqname);

			// ���� ���� ����: ���̿� Whole Sequence �Է� �ޱ�
			fscanf(fspace, "%d\n", &len);
			for (j = 0; j < MAX; j++)
				whole[j] = '\0';
			fscanf(fspace, "%s\n\n", whole);

			// Motif ���� : Motif ���� ������ �� ������ ���� �Է� �ޱ�
			fscanf(fmotif, "%d\n", &motif_type);
			motif = (char **)malloc(sizeof(char*) * motif_type);
			for (j = 0; j < motif_type; j++)
				motif[j] = (char *)malloc(sizeof(char) * (LENGTH + 1));
			for (j = 0; j < motif_type; j++) {
				for (k = 0; k <= LENGTH; k++)
					motif[j][k] = '\0';
			}
			for (j = 0; j < motif_type; j++)
				fscanf(fmotif, "%s %d\n", motif[j], &temp);		// ������ �Է¹ް� ������ temp�� ó���ؼ� ����

			// Interaction Ƚ�� ����ؼ� ���
			for (j = 0; j < motif_type; j++) {
				for (k = 0; k < motif_type; k++) {
					for (l = 0 ; l < LENGTH; l++)
						interaction[l] = motif[j][l];
					for (l = LENGTH; l < LENGTH * 2; l++)
						interaction[l] = motif[k][l - LENGTH];
					interaction[LENGTH * 2] = '\0';

					cnt = 0;
					for (l = 0; l <= len - (LENGTH * 2); l++) {
						flag = 0;	// 0�� ��� ����, 1�� ��� �ٸ�
						for (m = 0; m < LENGTH * 2; m++) {
							if ((interaction[m] != whole[l + m]) && (interaction[m] != '*')) {
								flag = 1;
								break;
							}
						}

						if (flag == 0)	// ���ٸ� ī��Ʈ
							cnt++;
					}

					if (cnt > 0) {
						fprintf(fint, "%s %d\n", interaction, cnt);
						total_int += cnt;
						total_type++;
					}
				}
			}

			fprintf(fsum, "\t>Total %d interaction(s) of %d interaction type(s)!\n", total_int, total_type);
			fprintf(fname, "%d ", total_type);	// Interaction DB ���� ��� ������Ʈ!

			for (j = 0; j < motif_type; j++)
				free(motif[j]);
			free(motif);
		}

		fprintf(fname, "\n");
		fclose(fmotif);
		fclose(fspace);
		fclose(fint);
	}

	fclose(fmotif_list);
	fclose(fspace_list);

	fclose(fsum);
	fclose(fname);
}