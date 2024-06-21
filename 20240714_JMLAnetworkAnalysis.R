#パッケージ準備
install.packages("igraph")
library(igraph)

#データ読み込み
g5 <- read.csv("20240213net5.csv")

#有向ネットワークグラフ作成@windows環境
g <- graph_from_data_frame(g5, directed=TRUE)
plot(g, layout=layout_with_fr, vertex.label.cex=0.8, 
     vertex.frame.color="white", edge.arrow.width=0.6, 
     vertex.size=4.5, edge.arrow.size=0.1, edge.curved=0.15)
#レイアウト固定。
coord <- layout_with_fr(g)

#デンドログラム作成
ceb <- cluster_edge_betweenness(g)
plot_dendrogram(ceb, mode="hclust")

#ネットワークグラフとデンドログラムの掛け合わせ
plot(ceb, g, layout=coord, vertex.label.cex=0.8, 
     vertex.frame.color="white", edge.arrow.width=0.6, vertex.size=4.5, 
     edge.arrow.size=0.1, edge.curved=0.15)

#推移性(global)とクラスター係数(average)計算
transitivity(g, type = "global")
transitivity(g, type = "average")

#次数カウント&ヒストグラム化&CSV出力
deg <- degree(g, mode="all")
deg_in <- degree(g, mode="in")
deg_out <- degree(g, mode="out")
hist(deg, breaks=14)
write.csv(deg, file="deg.csv")

#ネットワーク距離の計算とヒストグラム&CSV出力
network_distances <- distances(g)
hist(network_distances, breaks=10)
write.csv(network_distances, file="network_distances.csv")

#平均最短経路:無指向
mean_distance(g, directed=F)

#3点モチーフ数の計測
(count_motifs <- motifs(g, size=3))

#モチーフ一覧
#モチーフにカウントした数を書き込む@windows環境
co <- matrix(c(1,1, 0,0, 2,0), ncol =2, byrow = TRUE)
windows()
par(mfrow = c(4,4), mar = c(0,0,0,0))
for (i in 1:16) {
  plot(graph_from_isomorphism_class(3, i-1, directed = TRUE),
       layout = co, vertex.color = "red", vertex.frame.color = NA,
       vertex.label = NA, edge.color = "darkgray",
       edge.arrow.size = 1.5, edge.width = 4, frame = TRUE, margin = 0.1)
  text(0,0, count_motifs[i], cex =4, col = "blue")}
